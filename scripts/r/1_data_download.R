#팔요한 라이브러리 설치치
install.packages("BiocManager")
BiocManager::install("GEOquery")

#라이브러리 불러오기
library(GEOquery)

gsm_ids <- paste0("GSM", c(1868036:1868427, 1946528:1946557))  # GSM Sample ID(여러개 넣는거 가능)

download_dir <- "./" #이 디렉토리에 데이터 파일 저장(형식은 'soft')

for (gsm_id in gsm_ids) {
  tryCatch({
    gsm_data <- getGEO(gsm_id, destdir = download_dir) # GSM 한 샘플
    Sys.sleep(2) # 대기 시간 2초
    gsm_data_list[[gsm_id]] <- gsm_data # GSM 샘플 리스트
    print(paste0(gsm_id, "/", length(gsm_ids))) # 진행 과정 확인용용
  }, error = function(e) {
    print(paste0("Error: ", e$message))
  })
}
# 다운로드한 전체 샘플데이터 확인
View(gsm_data_list)