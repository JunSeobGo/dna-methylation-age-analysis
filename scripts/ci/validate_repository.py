"""원천 데이터 다운로드 없이 저장소 파일을 검증한다."""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROHIBITED_PREFIXES = ("data/raw/", "data/processed/", "outputs/")
PROHIBITED_SUFFIXES = (".zip", ".xlsx", ".hwp", ".pdf", ".pptx", ".ppt", ".docx", ".doc")
MARKDOWN_LINK_PATTERN = re.compile(r"\[[^\]]+\]\(([^)]+)\)")

CONFIG_SCHEMAS = {
    "cohort_registry.csv": {
        "required": {
            "series_id", "classification", "proposed_role", "expected_tissue",
            "expected_platform", "include_in_age_model", "reason", "source_url",
        },
        "key": ("series_id",),
    },
    "dataset_download_manifest.csv": {
        "required": {
            "dataset_id", "resource_id", "resource_type", "remote_url", "local_path",
            "expected_size_mb", "required_for_stage", "description",
        },
        "key": ("dataset_id", "resource_id"),
    },
    "external_validation_candidates.csv": {
        "required": {"accession", "priority", "status", "reason", "source_url"},
        "key": ("accession",),
    },
    "gse207605_cohort_manifest.csv": {
        "required": {
            "source_series_id", "expected_samples", "proposed_role", "model_policy",
            "reason", "remote_url", "source_url",
        },
        "key": ("source_series_id",),
    },
    "satsa_evaluation_policy.csv": {
        "required": {
            "policy_version", "dataset", "comparison", "older_age_min",
            "overall_mae_noninferiority_margin", "older_mae_improvement_min",
            "absolute_bias_worsening_margin", "bootstrap_replicates", "bootstrap_seed",
            "confidence_level", "cluster_unit",
        },
        "key": ("policy_version",),
    },
}


def tracked_files() -> list[str]:
    """Git이 추적하는 파일만 반환한다."""
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    return [path for path in result.stdout.decode("utf-8").split("\0") if path]


def validate_tracked_files() -> list[str]:
    """원천 데이터와 대용량 문서가 실수로 커밋되지 않았는지 검사한다."""
    violations: list[str] = []
    for path in tracked_files():
        normalized = path.replace("\\", "/")
        if normalized.startswith(PROHIBITED_PREFIXES):
            violations.append(f"추적 중인 원천 또는 생성 데이터: {path}")
        if normalized.lower().endswith(PROHIBITED_SUFFIXES):
            violations.append(f"추적 중인 대용량 문서: {path}")
    return violations


def validate_notebooks() -> list[str]:
    """노트북 파일이 유효한 JSON인지 검사한다."""
    violations: list[str] = []
    for notebook in (ROOT / "notebooks").glob("*.ipynb"):
        try:
            with notebook.open(encoding="utf-8") as file:
                json.load(file)
        except (OSError, json.JSONDecodeError) as error:
            violations.append(f"유효하지 않은 노트북 {notebook.relative_to(ROOT)}: {error}")
    return violations


def validate_markdown_links() -> list[str]:
    """Git이 추적하는 Markdown의 로컬 상대경로 링크를 검사한다."""
    violations: list[str] = []
    for relative_path in tracked_files():
        if not relative_path.lower().endswith(".md"):
            continue
        document = ROOT / relative_path
        try:
            content = document.read_text(encoding="utf-8")
        except OSError as error:
            violations.append(f"Markdown 파일을 읽을 수 없음 {relative_path}: {error}")
            continue

        for raw_target in MARKDOWN_LINK_PATTERN.findall(content):
            target = raw_target.strip().strip("<>")
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            path_only = target.split("#", maxsplit=1)[0]
            if not path_only:
                continue
            linked_path = document.parent / path_only
            if not linked_path.exists():
                violations.append(f"깨진 Markdown 링크 {relative_path}: {target}")
    return violations


def read_csv(relative_path: Path) -> tuple[list[str], list[dict[str, str]], list[str]]:
    """CSV를 읽고 헤더, 행, 구조 위반을 반환한다."""
    violations: list[str] = []
    path = ROOT / relative_path
    try:
        with path.open(encoding="utf-8-sig", newline="") as file:
            reader = csv.DictReader(file)
            headers = reader.fieldnames or []
            if not headers:
                return [], [], [f"헤더가 없는 CSV: {relative_path}"]
            if len(headers) != len(set(headers)):
                violations.append(f"중복 헤더가 있는 CSV: {relative_path}")
            if any(not header.strip() for header in headers):
                violations.append(f"빈 헤더가 있는 CSV: {relative_path}")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        return [], [], [f"CSV를 읽을 수 없음 {relative_path}: {error}"]

    if not rows:
        violations.append(f"데이터 행이 없는 CSV: {relative_path}")
    for number, row in enumerate(rows, start=2):
        if None in row:
            violations.append(f"헤더보다 값이 많은 CSV 행 {relative_path}:{number}")
        if all(not (value or "").strip() for key, value in row.items() if key is not None):
            violations.append(f"완전히 빈 CSV 행 {relative_path}:{number}")
    return headers, rows, violations


def validate_unique_keys(
    relative_path: Path,
    rows: list[dict[str, str]],
    keys: tuple[str, ...],
) -> list[str]:
    """설정 CSV의 식별 키가 비어 있거나 중복되는지 검사한다."""
    violations: list[str] = []
    seen: set[tuple[str, ...]] = set()
    for number, row in enumerate(rows, start=2):
        value = tuple((row.get(key) or "").strip() for key in keys)
        if any(not item for item in value):
            violations.append(f"빈 식별 키 {relative_path}:{number} ({', '.join(keys)})")
        elif value in seen:
            violations.append(f"중복 식별 키 {relative_path}:{number} ({', '.join(value)})")
        seen.add(value)
    return violations


def validate_config_values(name: str, rows: list[dict[str, str]]) -> list[str]:
    """분석 정책에 쓰이는 핵심 값의 허용 범위와 형식을 검사한다."""
    violations: list[str] = []
    allowed_values = {
        "cohort_registry.csv": {
            "classification": {"baseline", "candidate", "exclude"},
            "include_in_age_model": {"yes", "no", "conditional"},
        },
        "gse207605_cohort_manifest.csv": {
            "model_policy": {"include", "exclude", "conditional"},
        },
        "dataset_download_manifest.csv": {
            "required_for_stage": {"yes", "no", "conditional", "external_validation"},
        },
    }
    for number, row in enumerate(rows, start=2):
        for column, allowed in allowed_values.get(name, {}).items():
            value = (row.get(column) or "").strip()
            if value not in allowed:
                violations.append(f"허용되지 않은 값 {name}:{number} {column}={value!r}")

        for column in ("remote_url", "source_url"):
            if column in row and not (row.get(column) or "").startswith("https://"):
                violations.append(f"HTTPS URL이 아님 {name}:{number} {column}")

        if name == "dataset_download_manifest.csv":
            local_path = (row.get("local_path") or "").replace("\\", "/")
            if not local_path.startswith("data/raw/") or ".." in Path(local_path).parts:
                violations.append(f"허용되지 않은 원천 데이터 경로 {name}:{number}")
            try:
                if float(row.get("expected_size_mb") or "") <= 0:
                    raise ValueError
            except ValueError:
                violations.append(f"양수가 아닌 expected_size_mb {name}:{number}")
    return violations


def validate_satsa_policy(rows: list[dict[str, str]]) -> list[str]:
    """외부 평가 전에 고정한 SATSA 판정 기준이 우발적으로 바뀌지 않았는지 검사한다."""
    if len(rows) != 1:
        return ["satsa_evaluation_policy.csv는 정확히 한 행이어야 함"]

    expected = {
        "policy_version": "satsa-age-density-v1",
        "dataset": "S-BSST1206",
        "older_age_min": "80",
        "overall_mae_noninferiority_margin": "0.20",
        "older_mae_improvement_min": "0.50",
        "absolute_bias_worsening_margin": "0.50",
        "bootstrap_replicates": "10000",
        "bootstrap_seed": "20260804",
        "confidence_level": "0.95",
        "cluster_unit": "TID",
    }
    return [
        f"고정된 SATSA 정책 값이 다름: {column}={rows[0].get(column)!r} (기대값 {value!r})"
        for column, value in expected.items()
        if (rows[0].get(column) or "").strip() != value
    ]


def validate_config_csvs() -> list[str]:
    """config 폴더의 모든 CSV와 핵심 스키마를 원천 데이터 없이 검사한다."""
    violations: list[str] = []
    config_dir = ROOT / "config"
    for path in sorted(config_dir.glob("*.csv")):
        relative_path = path.relative_to(ROOT)
        headers, rows, csv_violations = read_csv(relative_path)
        violations.extend(csv_violations)

        schema = CONFIG_SCHEMAS.get(path.name)
        if not schema or not headers:
            continue
        missing = sorted(schema["required"] - set(headers))
        if missing:
            violations.append(f"필수 열 누락 {relative_path}: {', '.join(missing)}")
            continue
        violations.extend(validate_unique_keys(relative_path, rows, schema["key"]))
        violations.extend(validate_config_values(path.name, rows))
        if path.name == "satsa_evaluation_policy.csv":
            violations.extend(validate_satsa_policy(rows))
    return violations


def main() -> int:
    violations = (
        validate_tracked_files()
        + validate_notebooks()
        + validate_markdown_links()
        + validate_config_csvs()
    )
    if violations:
        print("저장소 검증에 실패했습니다:", file=sys.stderr)
        print("\n".join(f"- {violation}" for violation in violations), file=sys.stderr)
        return 1

    print("저장소 검증을 통과했습니다.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
