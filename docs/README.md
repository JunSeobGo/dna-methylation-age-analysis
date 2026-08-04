# 프로젝트 문서 안내

이 폴더는 분석 절차, 결과 해석, 기술·분석 의사결정을 관리합니다. 처음 보는 사람은 아래 순서로 읽으면 현재 상태와 재현 방법을 가장 빠르게 파악할 수 있습니다.

## 먼저 읽을 문서

1. [프로젝트 개요](PROJECT_OVERVIEW.md): 목표, 현재 모델 상태, 전체 분석 흐름
2. [다중 코호트 학습 풀](multicohort-training-pool.md): 포함·조건부·제외 기준과 5,541명 학습 bundle
3. [최종 모델과 외부 검증](final-model-and-external-validation.md): 표준 모델 잠금과 GSE87571 결과
4. [연령 편향 완화 실험](age-bias-mitigation-experiment.md): age-density 후보 선정 근거
5. [후보 모델 잠금과 새 외부 검증](age-density-model-lock-and-next-validation.md): 현재 후보의 적용 경계와 SATSA 계획

## 문서 분류

| 분류 | 문서 | 역할 |
| --- | --- | --- |
| 데이터 선정 | [외부 코호트 선정](cohort-selection.md), [GSE40279 준비](gse40279-data-preparation.md), [다중 코호트 학습 풀](multicohort-training-pool.md) | 코호트 포함·제외와 데이터 준비 근거 |
| 모델 검증 | [최종 모델과 외부 검증](final-model-and-external-validation.md), [학습곡선](learning-curve-analysis.md) | 일반화 성능과 표본 수 병목 점검 |
| 편향 진단 | [연령·코호트 편향 감사](age-cohort-bias-audit.md), [연령 편향 완화 실험](age-bias-mitigation-experiment.md) | 취약 연령대 진단과 후보 비교 |
| 현재 계획 | [후보 모델 잠금과 새 외부 검증](age-density-model-lock-and-next-validation.md) | 외부 검증 전 잠금 상태와 다음 단계 |
| 평가 정책 | [SATSA 잠금 모델 비교 평가](satsa-locked-evaluation.md) | 새 외부 결과 확인 전 고정한 지표·교체 조건·가족 bootstrap |
| 결과 보고 | [SATSA 외부 검증 보고서 템플릿](satsa-evaluation-report-template.md) | 실제 결과를 선택적으로 누락하지 않도록 사전에 고정한 보고 항목 |
| 재현 환경 | [R 실행환경 재현 기준](r-environment-reproducibility.md) | 잠금 모델과 일치해야 하는 R·glmnet 버전과 점검 방법 |
| 저장소 운영 | [GitHub 운영 가이드](github-workflow.md), [기여 가이드](../CONTRIBUTING.md) | Issue·브랜치·PR·커밋 규칙 |
| 의사결정 | [ADR 인덱스](adr/README.md) | 변경하기 어려운 기술·분석 결정과 이유 |

## 코드 분류

- `scripts/r/<번호>_*`: 현재 재현 파이프라인. 번호 순서가 데이터 준비부터 후보 잠금과 외부 평가까지의 실행 순서입니다.
- `scripts/r/Methyl_*`, `scripts/r/GSM샘플 저장.R`: 초기 분석과 탐색 과정의 레거시 스크립트입니다. 현행 모델 성능 재현에는 사용하지 않습니다.
- `notebooks/`: 초기 Python 전처리·탐색 분석 기록입니다. 현행 Elastic Net 평가 파이프라인과 구분합니다.
- `scripts/ci/`: 원천 데이터나 대용량 파일이 Git에 포함되지 않았는지 검사합니다.

## 문서 작성 규칙

- 실행 방법과 분석 결과는 해당 주제 문서에 기록합니다.
- 장기간 유지할 기술·분석 선택은 `docs/adr/`에 ADR로 기록합니다.
- 수치가 바뀔 수 있는 결과는 ADR에 복제하지 않고 결과 문서를 링크합니다.
- 개인 면접 답변, 회고 초안, 비공개 메모는 `portfolio-private/`에 두며 GitHub에 올리지 않습니다.
- 원천 데이터와 재생성 가능한 결과는 `data/`, `outputs/`에 두고 Git에서 제외합니다.
