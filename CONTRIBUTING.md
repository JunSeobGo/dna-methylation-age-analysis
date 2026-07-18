# 기여 및 작업 가이드

개인 프로젝트이지만 분석의 재현성과 변경 이력의 검토 가능성을 위해 가벼운 협업 흐름을 사용합니다.

## 브랜치 전략

| 브랜치 | 역할 |
| --- | --- |
| `main` | 문서화와 검증이 끝난 안정 버전 |
| `develop` | 다음 안정 버전을 위한 통합 브랜치 |
| `feature/<name>` | 하나의 기능 또는 분석 개선 작업 |
| `fix/<name>` | 재현성, 계산, 데이터 처리 오류 수정 |
| `docs/<name>` | 문서만 변경하는 작업 |

`main`과 `develop`에는 직접 커밋하지 않습니다. 모든 작업은 `develop`에서 분기한 작업 브랜치에서 수행하고 Pull Request로 `develop`에 병합합니다. 안정화 시에는 `develop`에서 `main`으로 별도 Pull Request를 생성합니다.

```bash
git switch develop
git pull --ff-only
git switch -c feature/<작업-이름>
```

## 작업 흐름

1. GitHub Issue에 목표, 작업 범위, 완료 조건을 작성합니다.
2. `develop`에서 하나의 작업 브랜치를 생성합니다.
3. 작은 단위로 커밋하고 검증합니다.
4. `develop`을 대상으로 Pull Request를 열고 검증 결과를 기록합니다.
5. 변경 파일과 GitHub Actions 결과를 확인한 뒤 병합합니다.

## 커밋 메시지

Conventional Commits 형식을 사용합니다.

```text
feat: 코호트 메타데이터 검증 추가
fix: 누락된 CpG probe 처리 수정
docs: 데이터 다운로드 절차 문서화
refactor: clock 계산 모듈 분리
test: 노트북 JSON 검증 추가
chore: 로컬 개발 가이드 정리
ci: 저장소 검증 workflow 추가
```

## Pull Request 전 점검

```bash
python scripts/ci/validate_repository.py
git diff --check
```

원천·처리 데이터, 생성 결과, API 키, 개인식별 가능 정보는 커밋하지 않습니다. `data/raw/`, `data/processed/`, `outputs/`는 로컬에서만 관리합니다.
