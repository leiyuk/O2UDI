from __future__ import annotations

from extract_cosmos_debris_norad_ids import extract_cosmos_debris_norad_ids, format_python_list
from pathlib import Path

_, cosmos_2251_20100101 = extract_cosmos_debris_norad_ids(Path(r"tle_snapshot_output/tle_snapshot_20090320_000000_UTC.3le"))
_, cosmos_2251_20200101 = extract_cosmos_debris_norad_ids(Path(r"tle_snapshot_output/tle_snapshot_20190320_000000_UTC.3le"))
print(len(cosmos_2251_20100101), len(cosmos_2251_20200101))
cosmos_2251_20100101_20200101 = list(set(cosmos_2251_20100101) & set(cosmos_2251_20200101))
print(len(cosmos_2251_20100101_20200101))
# print(sorted(cosmos_2251_20100101_20200101))
