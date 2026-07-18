# R Scripts

This folder contains the R analysis scripts used for the project.

## 다중 코호트 학습 풀 감사

`6_audit_gse207605_training_pool.R`은 `GSE207605`의 19개 구성 코호트를 내려받고 표본 수, 연령 범위, beta 값 범위, GSM 중복, 1차 학습 후보의 공통 CpG를 검사한다.

```powershell
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/6_audit_gse207605_training_pool.R --dry-run
& "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" scripts/r/6_audit_gse207605_training_pool.R --download
```

코호트별 포함·조건부·제외 판단은 `config/gse207605_cohort_manifest.csv`에서 관리하며, 상세 원칙은 `docs/multicohort-training-pool.md`를 참고한다.

## Files
- 1_data_download.R: download GEO sample metadata
- 1_data_download_from_geo.r: duplicate download script retained for convenience
- GSM샘플 저장.R: sample metadata saving helper
- Methyl_age_230526.R: DNA methylation age modeling workflow
- Methyl_시각화(Ridge, Elastic net, Lasso).R: visualization and modeling experiments
