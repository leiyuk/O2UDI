#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从指定的 .tle / .3le 文件中提取以下两类碎片的 NORAD ID：

1. COSMOS 1408 DEB
2. COSMOS 2251 DEB

适用文件格式：
    0 COSMOS 1408 DEB
    1 49516U ...
    2 49516  ...

或名称行不带前导 0：
    COSMOS 1408 DEB
    1 49516U ...
    2 49516  ...

说明：
- 仅提取名称中明确包含 “DEB” 的碎片，不会包含原始卫星。
- 提取结果自动去重并按 NORAD ID 从小到大排序。
- 如果文件只有两行根数、没有目标名称行，则无法仅凭 TLE 判断碎片来源。
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

# =============================================================================
# 全局配置
# =============================================================================

# 是否将提取结果额外保存为一个 Python 文件
SAVE_RESULT_TO_PY = True

# 结果文件路径
OUTPUT_PY_FILE = Path("cosmos_debris_norad_ids.py")

# =============================================================================
# 内部常量
# =============================================================================

COSMOS_1408_NAME = "COSMOS 1408 DEB"
COSMOS_2251_NAME = "COSMOS 2251 DEB"

# 匹配 TLE 第一行开头，例如：
# 1 49516U  ...
# 1 33757U  ...
TLE_LINE1_PATTERN = re.compile(r"^1\s+(\d{1,6})")


def read_text_with_fallback(file_path: Path) -> str:
    """使用常见编码读取文本文件。"""
    encodings = ("utf-8-sig", "utf-8", "gbk", "latin-1")
    last_error: Optional[UnicodeDecodeError] = None

    for encoding in encodings:
        try:
            return file_path.read_text(encoding=encoding)
        except UnicodeDecodeError as exc:
            last_error = exc

    raise RuntimeError(f"无法识别文件编码：{file_path}") from last_error


def normalize_object_name(line: str) -> str:
    """
    规范化 TLE 名称行。

    例如：
        "0 COSMOS 1408 DEB" -> "COSMOS 1408 DEB"
        "COSMOS   1408 DEB" -> "COSMOS 1408 DEB"
    """
    name = line.strip()

    if name.startswith("0 "):
        name = name[2:].strip()

    # 合并连续空白，并统一为大写
    return " ".join(name.upper().split())


def classify_debris_name(name: str) -> Optional[str]:
    """
    根据名称判断碎片所属类别。

    返回：
        "cosmos_1408"
        "cosmos_2251"
        None
    """
    normalized = normalize_object_name(name)

    if COSMOS_1408_NAME in normalized:
        return "cosmos_1408"

    if COSMOS_2251_NAME in normalized:
        return "cosmos_2251"

    return None


def extract_norad_id_from_line1(line1: str) -> int:
    """从 TLE 第一行提取 NORAD_CAT_ID。"""
    match = TLE_LINE1_PATTERN.match(line1.strip())

    if match is None:
        raise ValueError(f"无法从 TLE 第一行提取 NORAD ID：{line1!r}")

    return int(match.group(1))


def extract_cosmos_debris_norad_ids(tle_file: Path, ) -> tuple[list[int], list[int]]:
    """
    从 TLE/3LE 文件中提取 COSMOS 1408 和 COSMOS 2251 碎片的 NORAD ID。

    文件应包含名称行。函数会把名称行与其后的 TLE 第一行对应起来。
    """
    if not tle_file.exists():
        raise FileNotFoundError(f"输入文件不存在：{tle_file.resolve()}")

    text = read_text_with_fallback(tle_file)
    lines = [line.strip() for line in text.splitlines() if line.strip()]

    cosmos_1408_ids: set[int] = set()
    cosmos_2251_ids: set[int] = set()

    current_object_name: Optional[str] = None
    name_line_count = 0
    tle_line1_count = 0
    unassociated_line1_count = 0

    for line in lines:
        # TLE 第一行
        if line.startswith("1 "):
            tle_line1_count += 1

            if current_object_name is None:
                unassociated_line1_count += 1
                continue

            debris_type = classify_debris_name(current_object_name)

            if debris_type is not None:
                norad_id = extract_norad_id_from_line1(line)

                if debris_type == "cosmos_1408":
                    cosmos_1408_ids.add(norad_id)
                elif debris_type == "cosmos_2251":
                    cosmos_2251_ids.add(norad_id)

            # 当前名称只对应紧随其后的这一组 TLE
            current_object_name = None
            continue

        # TLE 第二行本身不参与识别
        if line.startswith("2 "):
            continue

        # 其余非空行视为目标名称行
        current_object_name = line
        name_line_count += 1

    if name_line_count == 0:
        raise RuntimeError("文件中没有发现目标名称行。" "只有两行根数的 2LE 文件无法仅凭 TLE 判断碎片属于哪次解体事件。")

    cosmos_1408_list = sorted(cosmos_1408_ids)
    cosmos_2251_list = sorted(cosmos_2251_ids)

    print(f"[input] 文件：{tle_file.resolve()}")
    print(f"[input] 名称行数量：{name_line_count}")
    print(f"[input] TLE 第一行数量：{tle_line1_count}")

    if unassociated_line1_count:
        print(f"[warning] 有 {unassociated_line1_count} 条 TLE 第一行" "未找到对应的名称行，已跳过。")

    return cosmos_1408_list, cosmos_2251_list


def format_python_list(variable_name: str, values: list[int]) -> str:
    """将整数列表格式化为易读的 Python 代码。"""
    if not values:
        return f"{variable_name} = []"

    lines = [f"{variable_name} = ["]
    current_line = "    "

    for value in values:
        item = f"{value}, "

        if len(current_line) + len(item) > 100:
            lines.append(current_line.rstrip())
            current_line = "    " + item
        else:
            current_line += item

    if current_line.strip():
        lines.append(current_line.rstrip())

    lines.append("]")
    return "\n".join(lines)


def save_python_result(
    output_file: Path,
    cosmos_1408_ids: list[int],
    cosmos_2251_ids: list[int],
) -> None:
    """将两个列表保存为可直接导入的 Python 文件。"""
    output_file.parent.mkdir(parents=True, exist_ok=True)

    content = ("# -*- coding: utf-8 -*-\n"
               '"""从 TLE 文件中提取的 COSMOS 碎片 NORAD ID。"""\n\n' + format_python_list(
                   "cosmos_1408_debris_norad_ids",
                   cosmos_1408_ids,
               ) + "\n\n" + format_python_list(
                   "cosmos_2251_debris_norad_ids",
                   cosmos_2251_ids,
               ) + "\n")

    output_file.write_text(content, encoding="utf-8")


def main(input_file: Path) -> None:
    cosmos_1408_ids, cosmos_2251_ids = (extract_cosmos_debris_norad_ids(input_file))

    print()
    # print(format_python_list(
    #     "cosmos_1408_debris_norad_ids",
    #     cosmos_1408_ids,
    # ))
    print(f"\nCOSMOS 1408 碎片数量：{len(cosmos_1408_ids)}")

    print()
    # print(format_python_list(
    #     "cosmos_2251_debris_norad_ids",
    #     cosmos_2251_ids,
    # ))
    print(f"\nCOSMOS 2251 碎片数量：{len(cosmos_2251_ids)}")

    if SAVE_RESULT_TO_PY:
        save_python_result(
            OUTPUT_PY_FILE,
            cosmos_1408_ids,
            cosmos_2251_ids,
        )
        print(f"\n[output] Python 结果文件：{OUTPUT_PY_FILE.resolve()}")


if __name__ == "__main__":
    # 输入 TLE 文件路径
    INPUT_TLE_FILE = Path(r"tle_snapshot_output/tle_snapshot_20220301_000000_UTC.3le")
    main(input_file=INPUT_TLE_FILE)
