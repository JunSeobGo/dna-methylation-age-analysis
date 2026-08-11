# 필요한 패키지 로드
library(readr)
library(progress)

# 파일 경로 지정
file_path <- "C:/Users/ko911/downloads/GSE72680_beta_values.txt"

# 진행률 바 생성
pb <- progress_bar$new(total = 2)  # 파일 읽기와 View() 두 단계로 진행되므로 total은 2입니다

# 파일 읽기 및 실시간으로 진행률 업데이트
data <- read.table(file_path, header = TRUE, sep = "\t", progress = pb$tick())

# View() 함수 호출
View(data)

data<-read.table("C:/Users/ko911/downloads/GSE72680_beta_values.txt", header = TRUE, sep = "\t") %>% View()



#-------------------------------------------------------------------------------
# GSM 샘플 데이터 불러오기
library(GEOquery)

#gsm_ids <- paste0("GSM",1868036:1946557)  
gsm_ids <- paste0("GSM", c(1868036:1868427, 1946528:1946557))  # GSM ID(여러개 넣는거 가능)

download_dir <- "C:/Users/ko911/downloads" #이 디렉토리에 데이터 파일 저장(형식은 'soft')

gsm_data_list <- list()
gsm_ch <- list()

for (gsm_id in gsm_ids) {
  gsm_data <- getGEO(gsm_id, destdir = download_dir) # GSM 한 샘플
  gsm_data_list[[gsm_id]] <- gsm_data # # GSM 샘플 리스트
  gsm_ch[[gsm_id]] <- gsm_data@header[["characteristics_ch1"]] # 필요한 부분 추출(샘플 특징)
  print(paste0(gsm_id, "/", length(gsm_ids)))
} # 현재 1868473~1868499랑 1869310~ 부분 저장해야함.

# 다운로드한 전체 샘플데이터 확인
View(gsm_data_list)

# 데이터 특징 불러오기
gsm_ch # 
gsm_ch[["GSM1869310"]]
