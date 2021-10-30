# initialization -------------
library(googledrive)

f.google.drive.model.retrieve <- function(){
  
  # # model_030621.rds
  # drive_download(file="https://drive.google.com/file/d/1FWAm4xzzAWdWXymPAfTPpkVBsPci6ytK/view?usp=sharing",
  #                path="manuscript/models/model_030621.rds")
  
  # # model_060221.rds
  # drive_download(file="https://drive.google.com/file/d/1y5MxKfz1XAtp6nlvnOtjzAaKQt9BHssA/view?usp=sharing",
  #                path="manuscript/models/model_060221.rds")
  
  # "model_070621.rds"
  drive_download(file="https://drive.google.com/file/d/1jKaW6LYIR4J7k7grZRF9vYPjrzekFtHG/view?usp=sharing",
                 path="manuscript/models/model_070621.rds")
  
}  

f.google.drive.model.retrieve()

# end ----------