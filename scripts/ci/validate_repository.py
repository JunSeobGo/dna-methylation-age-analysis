"""원천 데이터 다운로드 없이 저장소 파일을 검증한다."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PROHIBITED_PREFIXES = ("data/raw/", "data/processed/", "outputs/")
PROHIBITED_SUFFIXES = (".zip", ".xlsx", ".hwp", ".pdf", ".pptx", ".ppt", ".docx", ".doc")


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


def main() -> int:
    violations = validate_tracked_files() + validate_notebooks()
    if violations:
        print("저장소 검증에 실패했습니다:", file=sys.stderr)
        print("\n".join(f"- {violation}" for violation in violations), file=sys.stderr)
        return 1

    print("저장소 검증을 통과했습니다.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
