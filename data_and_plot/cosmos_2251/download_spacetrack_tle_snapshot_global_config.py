#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 Space-Track 获取指定目标时刻的历史 TLE 快照。

筛选规则: 
1. 查询 [目标时刻 - lookback_days, 目标时刻] 范围内的 GP_HISTORY。
2. 对每个 NORAD_CAT_ID，仅保留满足 EPOCH <= 目标时刻的最新一条 TLE。
3. 本脚本不使用 SGP4，也不会把 TLE 推进到目标时刻。
4. 输出: 
   - CSV: 包含目标名称、NORAD 编号、TLE 历元、TLE 年龄、两行根数；
   - 3LE 文本: 名称行 + TLE line 1 + TLE line 2。

依赖: pip install requests

Windows PowerShell 设置账号密码: 
    $env:SPACETRACK_USERNAME="你的Space-Track账号"
    $env:SPACETRACK_PASSWORD="你的Space-Track密码"

使用方法:
1. 在脚本顶部“用户配置区”设置 TARGET_TIME_UTC、LOOKBACK_DAYS 和 OUTPUT_DIR。
2. 运行:
       python download_spacetrack_tle_snapshot_global_config.py
"""

from __future__ import annotations

import csv
import os
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator
from urllib.parse import quote

import requests

LOGIN_URL = "https://www.space-track.org/ajaxauth/login"
QUERY_ROOT = "https://www.space-track.org/basicspacedata/query"

# =============================================================================
# 用户配置区
# =============================================================================
# 目标 UTC 时刻。对每个 NORAD 目标，仅保留 EPOCH 不晚于该时刻的最新一条 TLE。
TARGET_TIME_UTC = "2019-05-10 00:00:00"

# 仅查询目标时刻之前这段天数内的历史 TLE。
# 若部分目标更新较慢，可适当增大，例如改为 14.0 或 30.0。
LOOKBACK_DAYS = 2.0

# 输出目录。
OUTPUT_DIR = "tle_snapshot_output"

# HTTP 连接超时和读取超时，单位为秒。
REQUEST_TIMEOUT = (30, 600)


def parse_utc_datetime(value: str) -> datetime:
    """解析常见 UTC 时间格式，返回带 UTC 时区的 datetime。"""
    text = value.strip().replace("Z", "+00:00")

    try:
        dt = datetime.fromisoformat(text)
    except ValueError as exc:
        raise ValueError(f"无法解析时间 {value!r}，建议格式: YYYY-MM-DD HH:MM:SS") from exc

    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)
    return dt


def format_api_time(dt: datetime) -> str:
    """转换为 Space-Track 查询路径使用的 UTC 时间字符串。"""
    return dt.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")


def build_query_url(target: datetime, lookback_days: float) -> str:
    """构造 GP_HISTORY CSV 查询地址。"""
    start = target - timedelta(days=lookback_days)
    epoch_range = f"{format_api_time(start)}--{format_api_time(target)}"
    encoded_range = quote(epoch_range, safe="-:")

    predicates = ",".join([
        "NORAD_CAT_ID",
        "OBJECT_NAME",
        "EPOCH",
        "TLE_LINE1",
        "TLE_LINE2",
    ])

    return (f"{QUERY_ROOT}/class/gp_history"
            f"/EPOCH/{encoded_range}"
            f"/orderby/NORAD_CAT_ID%20asc,EPOCH%20desc"
            f"/predicates/{predicates}"
            f"/format/csv"
            f"/emptyresult/show")


def nonempty_decoded_lines(response: requests.Response) -> Iterator[str]:
    """逐行解码 HTTP 响应，并去除空行与首行 BOM。"""
    first = True
    for raw_line in response.iter_lines(decode_unicode=True):
        if raw_line is None:
            continue

        if isinstance(raw_line, bytes):
            line = raw_line.decode("utf-8-sig", errors="replace")
        else:
            line = raw_line

        line = line.strip("\r\n")
        if not line:
            continue

        if first:
            line = line.lstrip("\ufeff")
            first = False
        yield line


def login(session: requests.Session, username: str, password: str) -> None:
    """登录 Space-Track，并将认证 Cookie 保存在 session 中。"""
    response = session.post(
        LOGIN_URL,
        data={
            "identity": username,
            "password": password
        },
        timeout=REQUEST_TIMEOUT,
    )
    response.raise_for_status()

    text = response.text.lower()
    if "login failed" in text or "invalid" in text and "password" in text:
        raise RuntimeError("Space-Track 登录失败，请检查账号和密码。")


def select_latest_tles_at_or_before_target(
    session: requests.Session,
    query_url: str,
    target: datetime,
) -> Dict[str, dict]:
    """
    流式读取结果，并为每个 NORAD_CAT_ID 保留 EPOCH <= 目标时刻的最新 TLE。
    """
    selected: Dict[str, dict] = {}
    rows_read = 0
    rows_skipped = 0

    with session.get(query_url, stream=True, timeout=REQUEST_TIMEOUT) as response:
        response.raise_for_status()

        content_type = response.headers.get("Content-Type", "").lower()
        if "text/html" in content_type:
            preview = response.text[:500]
            raise RuntimeError("查询返回了 HTML，而不是 CSV。可能是登录失效、查询地址错误或账号受限。\n" f"响应开头: {preview}")

        reader = csv.DictReader(nonempty_decoded_lines(response))
        if not reader.fieldnames:
            raise RuntimeError("Space-Track 返回了空响应。")

        required = {
            "NORAD_CAT_ID",
            "OBJECT_NAME",
            "EPOCH",
            "TLE_LINE1",
            "TLE_LINE2",
        }
        missing = required.difference(reader.fieldnames)
        if missing:
            raise RuntimeError(f"返回字段不完整，缺少: {sorted(missing)}；实际字段: {reader.fieldnames}")

        for row in reader:
            rows_read += 1
            norad_id = (row.get("NORAD_CAT_ID") or "").strip()
            line1 = (row.get("TLE_LINE1") or "").rstrip()
            line2 = (row.get("TLE_LINE2") or "").rstrip()
            epoch_text = (row.get("EPOCH") or "").strip()

            if not norad_id or not epoch_text or not line1 or not line2:
                rows_skipped += 1
                continue

            try:
                epoch = parse_utc_datetime(epoch_text)
            except ValueError:
                rows_skipped += 1
                continue

            if epoch > target:
                rows_skipped += 1
                continue

            previous = selected.get(norad_id)
            if previous is None or epoch > previous["_epoch_dt"]:
                row["_epoch_dt"] = epoch
                selected[norad_id] = row

    print(f"[query] 历史记录读取数: {rows_read}")
    print(f"[query] 跳过无效记录数: {rows_skipped}")
    print(f"[query] 最终目标数: {len(selected)}")
    return selected


def norad_sort_key(item: dict) -> tuple[int, str]:
    text = str(item.get("NORAD_CAT_ID", ""))
    try:
        return int(text), text
    except ValueError:
        return sys.maxsize, text


def write_outputs(
    records: Iterable[dict],
    target: datetime,
    output_dir: Path,
) -> tuple[Path, Path]:
    """输出 CSV 和 3LE 文本。"""
    output_dir.mkdir(parents=True, exist_ok=True)
    stamp = target.strftime("%Y%m%d_%H%M%S_UTC")
    csv_path = output_dir / f"tle_snapshot_{stamp}.csv"
    tle_path = output_dir / f"tle_snapshot_{stamp}.3le"

    ordered = sorted(records, key=norad_sort_key)

    # with csv_path.open("w", encoding="utf-8-sig", newline="") as csv_file, \
    #         tle_path.open("w", encoding="utf-8", newline="\n") as tle_file:
    with tle_path.open("w", encoding="utf-8", newline="\n") as tle_file:
        # fieldnames = [
        #     "NORAD_CAT_ID",
        #     "OBJECT_NAME",
        #     "EPOCH",
        #     "TLE_AGE_HOURS_AT_TARGET",
        #     "TLE_LINE1",
        #     "TLE_LINE2",
        # ]
        # writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        # writer.writeheader()

        for row in ordered:
            epoch = row["_epoch_dt"]
            age_hours = (target - epoch).total_seconds() / 3600.0
            name = (row.get("OBJECT_NAME") or f"NORAD {row['NORAD_CAT_ID']}").strip()
            line1 = row["TLE_LINE1"].rstrip()
            line2 = row["TLE_LINE2"].rstrip()

            # writer.writerow({
            #     "NORAD_CAT_ID": row["NORAD_CAT_ID"],
            #     "OBJECT_NAME": name,
            #     "EPOCH": epoch.isoformat().replace("+00:00", "Z"),
            #     "TLE_AGE_HOURS_AT_TARGET": f"{age_hours:.6f}",
            #     "TLE_LINE1": line1,
            #     "TLE_LINE2": line2,
            # })

            tle_file.write(f"0 {name}\n")
            tle_file.write(f"{line1}\n")
            tle_file.write(f"{line2}\n")

    return csv_path, tle_path


def main() -> None:
    # 账号密码仍通过环境变量读取，避免写入脚本。
    username = os.getenv("SPACETRACK_USERNAME", "").strip()
    password = os.getenv("SPACETRACK_PASSWORD", "").strip()
    if not username or not password:
        raise RuntimeError("未找到 Space-Track 账号密码。请设置环境变量 " "SPACETRACK_USERNAME 和 SPACETRACK_PASSWORD。")

    if LOOKBACK_DAYS <= 0:
        raise ValueError("全局变量 LOOKBACK_DAYS 必须大于 0。")

    target = parse_utc_datetime(TARGET_TIME_UTC)
    query_url = build_query_url(target, LOOKBACK_DAYS)

    print(f"[config] 目标时刻 UTC: {target.isoformat()}")
    print(f"[config] 回溯天数: {LOOKBACK_DAYS}")
    print(f"[config] 输出目录: {Path(OUTPUT_DIR).resolve()}")
    print(f"[config] 查询地址: {query_url}")

    session = requests.Session()
    session.headers.update({
        "User-Agent": "tle-snapshot-downloader/1.0",
        "Accept": "text/csv,application/json;q=0.9,*/*;q=0.8",
    })

    try:
        login(session, username, password)
        records = select_latest_tles_at_or_before_target(session, query_url, target)
    finally:
        session.close()

    if not records:
        raise RuntimeError("没有查询到任何 TLE。请检查目标时间，或适当增大全局变量 LOOKBACK_DAYS。")

    csv_path, tle_path = write_outputs(records.values(), target, Path(OUTPUT_DIR))

    ages = [(target - row["_epoch_dt"]).total_seconds() / 3600.0 for row in records.values()]
    print(f"[result] TLE 年龄最小值: {min(ages):.3f} h")
    print(f"[result] TLE 年龄最大值: {max(ages):.3f} h")
    print(f"[result] CSV: {csv_path.resolve()}")
    print(f"[result] 3LE: {tle_path.resolve()}")


if __name__ == "__main__":
    main()
