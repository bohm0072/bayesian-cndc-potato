# initialization -------------
library(googledrive)

if(file.exists("manuscript/models/model.rds")!=T){
  
  drive_download(file="https://drive.google.com/file/d/1FWAm4xzzAWdWXymPAfTPpkVBsPci6ytK/view?usp=sharing",
                 path="manuscript/models/model.rds")
  
}



# end ----------