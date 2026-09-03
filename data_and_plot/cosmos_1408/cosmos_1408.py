from __future__ import annotations

from extract_cosmos_debris_norad_ids import extract_cosmos_debris_norad_ids, format_python_list
from pathlib import Path

cosmos_1408_20220301, _ = extract_cosmos_debris_norad_ids(Path(r"tle_snapshot_output/tle_snapshot_20220515_000000_UTC.3le"))
cosmos_1408_20220601, _ = extract_cosmos_debris_norad_ids(Path(r"tle_snapshot_output/tle_snapshot_20220815_000000_UTC.3le"))
print(len(cosmos_1408_20220301), len(cosmos_1408_20220601))
cosmos_1408_20220301_20220601 = list(set(cosmos_1408_20220301) & set(cosmos_1408_20220601))
print(len(cosmos_1408_20220301_20220601))
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 Space-Track 获取多个 NORAD ID 在指定 UTC 时间段内的全部历史 TLE。

输出目录示例：
    tle_history_output/
        25544.txt
        33757.txt
        49516.txt

每个文件只保存对应 NORAD ID 的 TLE，按 EPOCH 升序排列。
仅当该目标去重后的有效 TLE 条目数严格大于 MIN_VALID_TLE_COUNT 时，
才会输出对应文件。

账号密码通过环境变量设置：

Windows PowerShell：
    $env:SPACETRACK_USERNAME="你的 Space-Track 账号"
    $env:SPACETRACK_PASSWORD="你的 Space-Track 密码"

依赖：
    pip install requests
"""

import csv
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Iterator
from urllib.parse import quote

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# =============================================================================
# 全局配置
# =============================================================================

# 指定多个 NORAD ID
NORAD_IDS = cosmos_1408_20220301_20220601.copy()

# 查询时间范围，按 UTC 解释
# 保留 START_TIME_UTC <= EPOCH <= END_TIME_UTC 的全部 TLE
START_TIME_UTC = "2022-05-15 00:00:00"
END_TIME_UTC = "2022-08-15 00:00:00"

# 输出目录；文件名为“NORAD ID.txt”
OUTPUT_DIR = Path("tle_history_output")

# 每批查询的 NORAD ID 数量
# 不对每颗卫星单独请求，而是用逗号分隔的 ID 批量查询
NORAD_BATCH_SIZE = 200

# 相邻两批查询之间的等待时间，单位为秒
REQUEST_INTERVAL_SECONDS = 5.0

# 有效 TLE 条目数筛选阈值。
# 只有去重后的有效 TLE 数量严格大于该值，才会输出该卫星的数据文件。
# 例如设置为 10，则至少需要 11 组有效 TLE 才会保留。
MIN_VALID_TLE_COUNT = 0

# 运行前是否删除本次 NORAD_IDS 对应的旧输出文件。
# 建议保持 True，避免上次运行中已输出、但本次未达到阈值的旧文件残留。
REMOVE_OLD_TARGET_FILES = True

# 相邻两组 TLE 之间是否增加空行
ADD_BLANK_LINE_BETWEEN_TLES = False

# 连接超时和读取超时，单位为秒
REQUEST_TIMEOUT = (30, 600)

# 429 或服务器临时错误时的最大重试次数
MAX_RETRIES = 3

# =============================================================================
# Space-Track 接口
# =============================================================================

LOGIN_URL = "https://www.space-track.org/ajaxauth/login"
QUERY_ROOT = "https://www.space-track.org/basicspacedata/query"


def parse_utc_datetime(value: str) -> datetime:
    """解析 UTC 时间字符串。"""
    text = value.strip().replace("Z", "+00:00")

    try:
        dt = datetime.fromisoformat(text)
    except ValueError as exc:
        raise ValueError(f"无法解析时间 {value!r}，建议格式：YYYY-MM-DD HH:MM:SS") from exc

    if dt.tzinfo is None:
        return dt.replace(tzinfo=timezone.utc)

    return dt.astimezone(timezone.utc)


def format_api_time(dt: datetime) -> str:
    """转换为 Space-Track 查询路径使用的 UTC 时间格式。"""
    return dt.astimezone(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")


def normalize_norad_ids(values: Iterable[int | str]) -> list[int]:
    """检查 NORAD ID，去重并按数值排序。"""
    result: set[int] = set()

    for value in values:
        try:
            norad_id = int(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"无效的 NORAD ID：{value!r}") from exc

        if norad_id <= 0:
            raise ValueError(f"NORAD ID 必须为正整数：{norad_id}")

        result.add(norad_id)

    if not result:
        raise ValueError("NORAD_IDS 不能为空。")

    return sorted(result)


def split_batches(values: list[int], batch_size: int) -> Iterator[list[int]]:
    """将 NORAD ID 拆分为多个批次。"""
    if batch_size <= 0:
        raise ValueError("NORAD_BATCH_SIZE 必须大于 0。")

    for start in range(0, len(values), batch_size):
        yield values[start:start + batch_size]


def build_query_url(
    norad_ids: list[int],
    start_time: datetime,
    end_time: datetime,
) -> str:
    """构造 GP_HISTORY 批量 CSV 查询地址。"""
    norad_text = ",".join(str(value) for value in norad_ids)
    encoded_norad_ids = quote(norad_text, safe=",")

    epoch_range = (f"{format_api_time(start_time)}--{format_api_time(end_time)}")
    encoded_epoch_range = quote(epoch_range, safe="-:")

    predicates = ",".join([
        "NORAD_CAT_ID",
        "OBJECT_NAME",
        "EPOCH",
        "TLE_LINE1",
        "TLE_LINE2",
    ])

    return (f"{QUERY_ROOT}/class/gp_history"
            f"/NORAD_CAT_ID/{encoded_norad_ids}"
            f"/EPOCH/{encoded_epoch_range}"
            f"/orderby/NORAD_CAT_ID%20asc,EPOCH%20asc"
            f"/predicates/{predicates}"
            f"/format/csv"
            f"/emptyresult/show")


def create_session() -> requests.Session:
    """创建带重试机制的 HTTP 会话。"""
    retry = Retry(
        total=MAX_RETRIES,
        connect=MAX_RETRIES,
        read=MAX_RETRIES,
        status=MAX_RETRIES,
        backoff_factor=2.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
        respect_retry_after_header=True,
        raise_on_status=False,
    )

    adapter = HTTPAdapter(max_retries=retry)

    session = requests.Session()
    session.mount("https://", adapter)
    session.headers.update({
        "User-Agent": "multi-norad-tle-history-downloader/1.0",
        "Accept": "text/csv,application/json;q=0.9,*/*;q=0.8",
    })
    return session


def login(
    session: requests.Session,
    username: str,
    password: str,
) -> None:
    """登录 Space-Track。"""
    response = session.post(
        LOGIN_URL,
        data={
            "identity": username,
            "password": password,
        },
        timeout=REQUEST_TIMEOUT,
    )
    response.raise_for_status()

    text = response.text.lower()
    if ("login failed" in text or ("invalid" in text and "password" in text)):
        raise RuntimeError("Space-Track 登录失败，请检查账号和密码。")


def nonempty_decoded_lines(response: requests.Response, ) -> Iterator[str]:
    """逐行解码响应，忽略空行并清除 BOM。"""
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


def download_batch(
    session: requests.Session,
    norad_ids: list[int],
    start_time: datetime,
    end_time: datetime,
) -> Dict[int, list[dict]]:
    """下载一批 NORAD ID 的历史 TLE。"""
    result: Dict[int, list[dict]] = {norad_id: [] for norad_id in norad_ids}
    requested_ids = set(norad_ids)
    query_url = build_query_url(norad_ids, start_time, end_time)

    print(f"[query] 目标数：{len(norad_ids)}，" f"范围：{norad_ids[0]} ... {norad_ids[-1]}")

    with session.get(
            query_url,
            stream=True,
            timeout=REQUEST_TIMEOUT,
    ) as response:
        response.raise_for_status()

        content_type = response.headers.get("Content-Type", "").lower()
        if "text/html" in content_type:
            preview = response.text[:500]
            raise RuntimeError("查询返回了 HTML，而不是 CSV。" "可能是登录失效、查询地址错误或账号受限。\n" f"响应开头：{preview}")

        reader = csv.DictReader(nonempty_decoded_lines(response))

        if not reader.fieldnames:
            return result

        required = {
            "NORAD_CAT_ID",
            "OBJECT_NAME",
            "EPOCH",
            "TLE_LINE1",
            "TLE_LINE2",
        }
        missing = required.difference(reader.fieldnames)

        if missing:
            fieldnames_text = ",".join(reader.fieldnames).upper()
            if "NO RESULTS" in fieldnames_text:
                return result

            raise RuntimeError(f"返回字段不完整，缺少：{sorted(missing)}；" f"实际字段：{reader.fieldnames}")

        for row in reader:
            norad_text = (row.get("NORAD_CAT_ID") or "").strip()
            epoch_text = (row.get("EPOCH") or "").strip()
            line1 = (row.get("TLE_LINE1") or "").rstrip()
            line2 = (row.get("TLE_LINE2") or "").rstrip()

            if not norad_text or not epoch_text or not line1 or not line2:
                continue

            try:
                norad_id = int(norad_text)
                epoch = parse_utc_datetime(epoch_text)
            except (TypeError, ValueError):
                continue

            if norad_id not in requested_ids:
                continue

            # 本地再次检查时间边界
            if epoch < start_time or epoch > end_time:
                continue

            result[norad_id].append({
                "OBJECT_NAME": ((row.get("OBJECT_NAME") or "").strip() or f"NORAD {norad_id}"),
                "EPOCH": epoch,
                "TLE_LINE1": line1,
                "TLE_LINE2": line2,
            })

    return result


def finalize_records(records: list[dict]) -> list[dict]:
    """按 EPOCH 和两行根数去重，并按 EPOCH 升序排序。"""
    unique: dict[tuple[datetime, str, str], dict] = {}

    for row in records:
        key = (
            row["EPOCH"],
            row["TLE_LINE1"],
            row["TLE_LINE2"],
        )
        unique[key] = row

    return sorted(
        unique.values(),
        key=lambda row: (
            row["EPOCH"],
            row["TLE_LINE1"],
            row["TLE_LINE2"],
        ),
    )


def write_norad_file(
    output_dir: Path,
    norad_id: int,
    records: list[dict],
) -> Path:
    """输出单个 NORAD ID 的 TLE 文件。"""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{norad_id}.txt"
    ordered = finalize_records(records)

    with output_path.open(
            "w",
            encoding="utf-8",
            newline="\n",
    ) as file:
        for index, row in enumerate(ordered):
            file.write(f"0 {row['OBJECT_NAME']}\n")
            file.write(f"{row['TLE_LINE1']}\n")
            file.write(f"{row['TLE_LINE2']}\n")

            if (ADD_BLANK_LINE_BETWEEN_TLES and index < len(ordered) - 1):
                file.write("\n")

    return output_path


def main() -> None:
    username = os.getenv("SPACETRACK_USERNAME", "").strip()
    password = os.getenv("SPACETRACK_PASSWORD", "").strip()

    if not username or not password:
        raise RuntimeError("未找到 Space-Track 账号密码。请设置环境变量 " "SPACETRACK_USERNAME 和 SPACETRACK_PASSWORD。")

    norad_ids = normalize_norad_ids(NORAD_IDS)
    start_time = parse_utc_datetime(START_TIME_UTC)
    end_time = parse_utc_datetime(END_TIME_UTC)

    if start_time > end_time:
        raise ValueError("START_TIME_UTC 不能晚于 END_TIME_UTC。")

    if REQUEST_INTERVAL_SECONDS < 0:
        raise ValueError("REQUEST_INTERVAL_SECONDS 不能小于 0。")

    if MIN_VALID_TLE_COUNT < 0:
        raise ValueError("MIN_VALID_TLE_COUNT 不能小于 0。")

    batches = list(split_batches(norad_ids, NORAD_BATCH_SIZE))
    all_records: Dict[int, list[dict]] = {norad_id: [] for norad_id in norad_ids}

    print(f"[config] NORAD ID 数量：{len(norad_ids)}")
    print(f"[config] UTC 起始时间：{start_time.isoformat()}")
    print(f"[config] UTC 结束时间：{end_time.isoformat()}")
    print(f"[config] 查询批次数：{len(batches)}")
    print(f"[config] TLE 条目筛选：有效条目数 > " f"{MIN_VALID_TLE_COUNT}")
    print(f"[config] 输出目录：{OUTPUT_DIR.resolve()}")

    session = create_session()

    try:
        login(session, username, password)
        print("[login] Space-Track 登录成功。")

        for index, batch_ids in enumerate(batches, start=1):
            print(f"[batch] {index}/{len(batches)}")

            batch_result = download_batch(
                session=session,
                norad_ids=batch_ids,
                start_time=start_time,
                end_time=end_time,
            )

            for norad_id, rows in batch_result.items():
                all_records[norad_id].extend(rows)

            if index < len(batches):
                print(f"[wait] 等待 {REQUEST_INTERVAL_SECONDS:.1f} 秒。")
                time.sleep(REQUEST_INTERVAL_SECONDS)

    finally:
        session.close()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 删除本次目标对应的旧文件，确保输出目录中不会残留
    # “本次有效条目数未超过阈值”的历史文件。
    if REMOVE_OLD_TARGET_FILES:
        removed_count = 0

        for norad_id in norad_ids:
            old_path = OUTPUT_DIR / f"{norad_id}.txt"

            if old_path.exists():
                old_path.unlink()
                removed_count += 1

        print(f"[cleanup] 已删除旧目标文件数：{removed_count}")

    output_tle_count = 0
    kept_ids: list[int] = []
    filtered_ids: list[tuple[int, int]] = []

    for norad_id in norad_ids:
        records = finalize_records(all_records[norad_id])
        valid_count = len(records)

        # 注意：这里是严格“大于”，不是大于等于。
        if valid_count > MIN_VALID_TLE_COUNT:
            kept_ids.append(norad_id)
            output_tle_count += valid_count

            output_path = write_norad_file(
                OUTPUT_DIR,
                norad_id,
                records,
            )
            # print(f"[output] {norad_id}：{valid_count} 组有效 TLE，" f"已保留 -> {output_path.resolve()}")
        else:
            filtered_ids.append((norad_id, valid_count))
            print(f"[filter] {norad_id}：{valid_count} 组有效 TLE，" f"未超过阈值 {MIN_VALID_TLE_COUNT}，不输出文件。")

    print()
    print(f"[result] 输入目标总数：{len(norad_ids)}")
    print(f"[result] 保留目标数：{len(kept_ids)}")
    print(f"[result] 筛除目标数：{len(filtered_ids)}")
    print(f"[result] 输出 TLE 总数：{output_tle_count}")
    print(f"[result] 输出目录：{OUTPUT_DIR.resolve()}")

    if kept_ids:
        # print(f"[result] 保留的 NORAD ID：{kept_ids}")
        pass

    if filtered_ids:
        print("[result] 被筛除的 NORAD ID 及有效条目数：" f"{filtered_ids}")


if __name__ == "__main__":
    main()
