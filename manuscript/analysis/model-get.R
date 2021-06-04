# initialization -------------
library(googledrive)



if(file.exists("manuscript/models/model.rds")!=T){
  
  # # model_030621.rds
  # drive_download(file="https://drive.google.com/file/d/1FWAm4xzzAWdWXymPAfTPpkVBsPci6ytK/view?usp=sharing",
  #                path="manuscript/models/model.rds")
  
  # model_060221.rds
  drive_download(file="https://drive.google.com/file/d/1y5MxKfz1XAtp6nlvnOtjzAaKQt9BHssA/view?usp=sharing",
                 path="manuscript/models/model.rds")
  
}

# end ----------