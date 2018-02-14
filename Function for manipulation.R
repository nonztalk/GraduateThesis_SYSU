ListFill <- function(list, max_length, std_list_content) {

  for (i in 1:length(list)) {
    if (length(list[[i]]) != max_length) {
      diff_pos <- match(setdiff(std_list_content, gsub(":.*", "", list[[i]])), std_list_content)
      list[[i]] <- append(list[[i]], values = paste(setdiff(std_list_content, gsub(":.*", "", list[[i]])), ": ", sep = ""))
    }
  }
  list <- lapply(list, StringSort, pattern = ":.*")
  return(list)
}

StringSort <- function(target, pattern) {
  require(stringr)
  sort.temp <- str_sort(gsub(pattern = pattern, "", target), locale = "en")
  target.finish <- target[order(match(gsub(pattern = pattern, "", target), sort.temp))]
  return(target.finish)
}

reshapeForBCData <- function(data, ch1, ch2) {
  
  if (! is.data.frame(data)) {
    stop("Not a data frame")
  }
  if (ch1 && ! ch2) {
    
    data.list <- strsplit(data$Ch1_characteristic, split = "1 = ")
    data.list <- lapply(data.list, StringSort, pattern = ":.*")
    std_content <- str_sort(unique(gsub(":.*", "", do.call(c, data.list))), locale = "en")
    data.list <- ListFill(data.list, max_length = length(std_content), std_list_content = std_content)
    
    chframe <- plyr::ldply(data.list)
    
    for (i in 1:ncol(chframe)) {
      colnames(chframe)[i] <- unique(substring(chframe[,i], 1, gregexpr(":", chframe[,i])[[1]][1] - 1))
      colnames(chframe)[i] <- gsub(" ", "_", colnames(chframe)[i])
      chframe[, i] <- gsub(";", "", chframe[, i])
      chframe[, i] <- gsub("^.*: ", "", chframe[, i])
    }
    
    chframe <<- chframe[! sapply(chframe, function(x) {all(x == "")})]
    
  }
  if (! ch1 && ch2) {
    
    data.list <- strsplit(data$Ch2_characteristic, split = "2 = ")
    data.list <- lapply(data.list, StringSort, pattern = ":.*")
    std_content <- str_sort(unique(gsub(":.*", "", do.call(c, data.list))), locale = "en")
    data.list <- ListFill(data.list, max_length = length(std_content), std_list_content = std_content)
    
    chframe <- plyr::ldply(data.list)
    
    for (i in 1:ncol(chframe)) {
      colnames(chframe)[i] <- unique(substring(chframe[,i], 1, gregexpr(":", chframe[,i])[[1]][1] - 1))
      colnames(chframe)[i] <- gsub(" ", "_", colnames(chframe)[i])
      chframe[, i] <- gsub(";", "", chframe[, i])
      chframe[, i] <- gsub("^.*: ", "", chframe[, i])
    }
    
    chframe <<- chframe[! sapply(chframe, function(x) {all(x == "")})]
    
  }
  
  if (ch1 && ch2) {
    
    data.list1 <- strsplit(data$Ch1_characteristic, split = "1 = ")
    data.list1 <- lapply(data.list1, StringSort, pattern = ":.*")
    std_content1 <- str_sort(unique(gsub(":.*", "", do.call(c, data.list1))), locale = "en")
    data.list1 <- ListFill(data.list1, max_length = length(std_content1), std_list_content = std_content1)
    
    data.list2 <- strsplit(data$Ch2_characteristic, split = "2 = ")
    data.list2 <- lapply(data.list2, StringSort, pattern = ":.*")
    std_content2 <- str_sort(unique(gsub(":.*", "", do.call(c, data.list2))), locale = "en")
    data.list2 <- ListFill(data.list2, max_length = length(std_content2), std_list_content = std_content2)
    
    chframe <- cbind(plyr::ldply(data.list1), plyr::ldply(data.list2))
    
    for (i in 1:ncol(chframe)) {
      colnames(chframe)[i] <- unique(substring(chframe[,i], 1, gregexpr(":", chframe[,i])[[1]][1] - 1))
      colnames(chframe)[i] <- gsub(" ", "_", colnames(chframe)[i])
      chframe[, i] <- gsub(";", "", chframe[, i])
      chframe[, i] <- gsub("^.*: ", "", chframe[, i])
    }
    
    chframe <<- chframe[! sapply(chframe, function(x) {all(x == "")})]
    
  }
  
}

ColMessage <- function(data, colname) {
  require(stringr)
  Col_Message <- c()
  for (col_name in colname) {
    if (length(unique(data[, col_name])) > 5){
      Col_content <- unique(data[, col_name])[1:5]
    }else{
      Col_content <- unique(data[, col_name])
    }
    Col_M <- str_c(col_name, ":", paste(Col_content, collapse = ","))
    Col_Message <- c(Col_Message, Col_M)
  }
  print(Col_Message)
}

readInput <- function() {
  cho <- readline(prompt = "Enter your choice:")
  cho <- strsplit(cho, split = " ")[[1]]
  return(as.integer(cho))
}

FillBCData <- function(target_data, target_col, target_col_position, fix = F) {
  require(stringr)
  filenames <- list.files(pattern = "^[0-9].*.csv")
  for (file in filenames) {
    from_data <- read.csv(file = file, header = T, stringsAsFactors = F)
    from_position <- grep(target_col, str_to_lower(colnames(from_data)[21:length(from_data)]), fixed = fix)
    if (length(from_position) == 0) {
      next
    }
    if (length(from_position) > 1) {
      ColMessage(from_data, colname = colnames(from_data)[20 + from_position])
      choice <- readInput()
      if (choice == -1) {
        next
      }
      if (choice == -2) {
        result_col <- apply(from_data[, colnames(from_data)[20 + from_position]], 1, paste, collapse = ";")
      }
      if (choice == -3) {
        result_col <- "Mark"
      }
      if (length(choice) > 1) {
        result_col <- apply(from_data[, colnames(from_data)[20 + from_position][choice]], 1, paste, collapse = ";")
      }
      result_col <- from_data[, gsub(":.*", "", colnames(from_data)[20 + from_position][choice])]
    }
    if (length(from_position) == 1) {
      result_col <- from_data[, 20 + from_position]
    }
    target_data[which(target_data$GSE_Summary == from_data$GSE_Summary), target_col_position] <- result_col
  }
  return(target_data)
}