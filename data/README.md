# Data Folder

This folder contains the project data files.

## Raw data
Place source data here before running the workflow.

## Processed data
Store intermediate or cleaned data outputs here.

## 현재 외부 검증 데이터 정책

- `GSE87571`은 기존 표준 모델의 잠금 외부 평가에 이미 사용했으므로 새 후보 선택에 재사용하지 않습니다.
- `data/raw/S-BSST1206/SATSA_pheno.txt`에는 공개 SATSA 표현형만 저장합니다.
- SATSA beta 약 6.7GB는 URL·해시·스트리밍 추출 절차를 확정하기 전에는 내려받지 않습니다.
- `data/processed/locked_age_density_candidate.rds`와 `satsa_external_sample_manifest.csv`는 로컬 재현 산출물이며 Git에서 제외합니다.

자세한 표본 규칙과 실행 순서는 `docs/age-density-model-lock-and-next-validation.md`를 참고하세요.
