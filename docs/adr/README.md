# Architecture Decision Records

ADR은 프로젝트의 중요한 기술·분석 결정을 당시 맥락과 함께 남기는 문서입니다. 결과 보고서가 “무엇이 나왔는지”를 설명한다면 ADR은 “왜 이 방식을 선택했고 무엇을 포기했는지”를 설명합니다.

## 상태

- `제안됨`: 논의 중이며 아직 파이프라인에 적용하지 않음
- `승인됨`: 현재 적용 중인 결정
- `대체됨`: 더 새로운 ADR이 대신함
- `폐기됨`: 더 이상 적용하지 않음

승인된 ADR의 결정을 바꿀 때는 기존 파일을 조용히 수정하지 않습니다. 새 ADR을 만들고 이전 ADR에 `대체됨: ADR-XXXX`를 표시합니다.

## 목록

| 번호 | 결정 | 상태 |
| --- | --- | --- |
| [ADR-0001](0001-version-control-boundary.md) | Git에는 코드·명세·문서만 저장 | 승인됨 |
| [ADR-0002](0002-cohort-grouped-validation.md) | 코호트 그룹 기반 중첩 교차검증 사용 | 승인됨 |
| [ADR-0003](0003-locked-external-validation.md) | 외부 검증은 모델 잠금 후 한 번 수행 | 승인됨 |
| [ADR-0004](0004-age-density-candidate.md) | age-density를 외부 검증 대기 후보로 잠금 | 승인됨 |
| [ADR-0005](0005-satsa-validation-design.md) | SATSA의 개인·가족 의존성과 대용량 beta 처리 정책 | 승인됨 |
| [ADR-0006](0006-satsa-model-replacement-criteria.md) | SATSA 모델 교체 조건을 결과 확인 전에 고정 | 승인됨 |
| [ADR-0007](0007-satsa-coverage-failure-response.md) | CpG coverage 미달로 869 트랙을 종결하고 610 제한 트랙 개시 | 승인됨 |

새 ADR은 [템플릿](0000-template.md)을 복사해 다음 번호로 작성합니다.
