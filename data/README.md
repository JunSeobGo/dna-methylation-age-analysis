# Data Folder

This folder contains the project data files.

## Raw data
Place source data here before running the workflow.

## Processed data
Store intermediate or cleaned data outputs here.

## 현재 외부 검증 데이터 정책

- `GSE87571`은 기존 표준 모델의 잠금 외부 평가에 이미 사용했으므로 새 후보 선택에 재사용하지 않습니다.
- `data/raw/S-BSST1206/SATSA_pheno.txt`에는 공개 SATSA 표현형만 저장합니다.
- SATSA beta 5개 파일의 URL·공식 바이트와 스트리밍 추출 절차는 고정했습니다. 내려받은 원본은 크기와 로컬 MD5를 검사합니다.
- `data/processed/locked_age_density_candidate.rds`와 `satsa_external_sample_manifest.csv`는 로컬 재현 산출물이며 Git에서 제외합니다.
- `satsa_external_validation_bundle.rds`는 사전 고정한 CpG coverage·결측·범위 게이트를 통과할 때만 생성합니다.

자세한 표본 규칙은 `docs/age-density-model-lock-and-next-validation.md`, 다운로드·추출 절차는 `docs/satsa-beta-extraction.md`를 참고하세요.
