library(shinythemes)
##require(XML)
library(xml2)
library(tidyverse)
library(plyr)
library(dplyr)
library(DT)
##library(Biostrings)

db_set <- c("Gmax_275_v2.0.fa",
            "Gmax_275_Wm82.a2.v1.protein_primaryTranscriptOnly.fa")
db_path <- c("~/database/Gmax_275_Wm82.a2.v1/assembly/")

ui <- fluidPage(theme = shinytheme("sandstone"),
                ##tagList(
                ##    tags$head(
                ##             tags$link(rel="stylesheet", type="text/css",href="style.css"),
                ##             tags$script(type="text/javascript", src = "busy.js")
                ##         )
                ##),
                sidebarLayout(
                    ##This block gives us all the inputs:
                    sidebarPanel(
                        #headerPanel('Shiny Blast'),
                        textAreaInput('query', 'Input sequence (FASTA):',
                                      value = "", placeholder = "",
                                      width = "100%", height="200px"),
                        selectInput("Reference", "Select Reference:",
                                    choices=c("William 82 a2v1", "W05"),
                                    width="300px"),
                        selectInput("dbtype", "Databse Type:",
                                    choices=c("Ref Type: Genome","Ref Type: Proteome"),
                                    selected = "Ref Type: Genome",
                                    width="300px"),
                        div(style="display:inline-block",
                            selectInput("program", "Program:",
                                        choices=c("blastn","tblastn"),
                                        width="148px")),
                        div(style="display:inline-block",
                            selectInput("eval", "e-value:",
                                        choices=c(1,0.001,1e-4,1e-5,1e-10),
                                        width="148px")),
                        textInput("word_size", "Word (W) length:",
                                  value = 11,
                                  placeholder = "default",
                                  width="148px"),
                        h5("Default = 11 for BLASTN, 3 for all others"),
                        textInput("num_alignments",
                                  "Max # of alignments to show:",
                                  value=50,
                                  placeholder="50",
                                  width="220px"),
                        fluidRow(
                            column(12,
                                   actionButton("blast", "GO", width = "120px"),
                                   align = "right"
                                   )),
                    width= 3),
                    
                    ##this snippet generates a progress indicator for long BLASTs
                    ##div(class = "busy",
                    ##    p("Calculation in progress.."), 
                    ##    img(src="https://i.stack.imgur.com/8puiO.gif", height = 100, width = 100,align = "center")
                    ##    ),
                
                    ##Basic results output
                    mainPanel(
                        h2("Blast Results"),
                        htmlOutput("program_version"),
                        htmlOutput("database_info"),
                        #br(),
                        hr(),
                        DT::dataTableOutput("blastResults"),
                        div(h3("Alignment:", style = "font-family: 'Lobster', cursive;"),
                            h4(htmlOutput("selected_alignment")),
                            verbatimTextOutput("alignment")),
                        width = 9)
                ),
                tags$head(tags$style(HTML("
                            #alignment {
                              font-size: 16px;
                              color: black;
                            }
                            #selected_alignment {
                               line-height: 1.5
                            }
                            ")))
)

extract_xml_value <- function(x, y){
    x %>% xml_find_all(paste0(".//",y))%>% xml_text()
}

Hsp_extract <- function(Hsp_xml){
    ## Hsp_xml is xml_nodeset of Hsp level
    ## Simple function to extract all value from Hsp level
    Hsp_num <- extract_xml_value(Hsp_xml, "Hsp_num")
    Hsp_score <- Hsp_xml %>% xml_find_all(".//Hsp_score") %>% xml_text()
    Hsp_bit_score <- Hsp_xml %>% xml_find_all(".//Hsp_bit-score") %>% xml_text()
    Hsp_evalue <- Hsp_xml %>% xml_find_all(".//Hsp_evalue") %>% xml_text()
    Hsp_align_len <- Hsp_xml %>% xml_find_all(".//Hsp_align-len") %>% xml_text()
    Hsp_identity <- Hsp_xml %>% xml_find_all(".//Hsp_identity") %>% xml_text()
    Hsp_gaps <- Hsp_xml %>% xml_find_all(".//Hsp_gaps") %>% xml_text()
    Hsp_positive <- Hsp_xml %>% xml_find_all(".//Hsp_positive") %>% xml_text()
    Hsp_query_from <- Hsp_xml %>% xml_find_all(".//Hsp_query-from") %>% xml_text()
    Hsp_query_to <- Hsp_xml %>% xml_find_all(".//Hsp_query-to") %>% xml_text()
    Hsp_hit_from <- Hsp_xml %>% xml_find_all(".//Hsp_hit-from") %>% xml_text()
    Hsp_hit_to <- Hsp_xml %>% xml_find_all(".//Hsp_hit-to") %>% xml_text()
    Hsp_query_frame <- Hsp_xml %>% xml_find_all(".//Hsp_query-frame") %>% xml_text()
    Hsp_hit_frame <- Hsp_xml %>% xml_find_all(".//Hsp_hit-frame") %>% xml_text()
    Hsp_qseq <- extract_xml_value(Hsp_xml, "Hsp_qseq")
    Hsp_hseq <- extract_xml_value(Hsp_xml, "Hsp_hseq")
    Hsp_midline <- extract_xml_value(Hsp_xml, "Hsp_midline")
    Hsp_record <- tibble(Hsp_num, Hsp_score, Hsp_bit_score, Hsp_evalue, Hsp_align_len, Hsp_identity, Hsp_gaps, Hsp_positive, Hsp_query_from, Hsp_query_to, Hsp_hit_from, Hsp_hit_to, Hsp_qseq, Hsp_hseq, Hsp_midline, Hsp_query_frame, Hsp_hit_frame)
    return(Hsp_record)
}
    
Parse_blastXML <- function(blastout){
    ## blastout is read_xml output xml2 object
    ##blastXML <- read_xml(blastXMLfile, as = "text")
    dataOut <- tibble()
    iter_query <- xml_find_all(blastout, "//Iteration")
    for(m in 1:length(iter_query)){
        query <- iter_query[m] %>%
            xml_find_all(".//Iteration_query-def") %>% xml_text()
        iter_hit <- iter_query[m] %>% xml_find_all(".//Hit")
        for(n in 1:length(iter_hit)){
            subject <- iter_hit[n] %>% xml_find_all(".//Hit_id") %>% xml_text()
            subject_len <- iter_hit[n] %>% extract_xml_value("Hit_len")
            iter_hsp <- iter_hit[n] %>% xml_find_all(".//Hsp")
            Hsp_list <- iter_hsp %>% map(Hsp_extract)
            dataOut <- dataOut %>%
                bind_rows(tibble(query, subject, subject_len, Hsp_list))
        }
    }
    return(dataOut)
}

Extract_blastXMLheader <- function(blastout){
    ## Get some program information
    ##blastout <- read_xml(blastXMLfile, as = "text")

    ## Use blastout (read_xml output) as input
    program_version <- blastout %>%
        extract_xml_value("BlastOutput_version") %>% unique()
    db_num <- blastout %>% 
        extract_xml_value("Statistics_db-num") %>% unique()
    db_len <- blastout %>%
        xml_find_all(".//Statistics_db-len") %>% xml_text() %>% unique()
    db_name <- blastout %>%
        xml_find_all(".//BlastOutput_db") %>% xml_text() %>% unique()
    queryID <- blastout %>%
        extract_xml_value("BlastOutput_query-def") %>% unique()
    queryLen <- blastout %>%
        extract_xml_value("BlastOutput_query-len") %>% unique()

    tibble(program_version, db_num, db_len, db_name, queryID, queryLen)
}

Alignment_singleLine <- function(Qstart, Qend, Sstart, Send, Qtext, Stext, middle){
    Qstart <- str_pad(Qstart, 15, "right")
    Qend <- str_pad(Qend, 15, "left")
    Sstart <- str_pad(Sstart, 15, "right")
    Send <- str_pad(Send, 15, "left")
    k1 <- str_pad("", 20, "both")
    k2 <- str_pad("", 15, "both")
    L1 <- paste0("Query  ", Qstart, Qtext, Qend, "\n")
    L2 <- paste0(k1, "  ", middle, k2, "\n")
    L3 <- paste0("Sbjct  ", Sstart, Stext, Send, "\n")
    paste0(L1, L2, L3)
}

Adjust_alignment <- function(Qstart, Qend, Sstart, Send, Qtext, Stext, middle, text_width){

    text_width <- as.numeric(text_width)
    start = 1
    end = start + text_width -1
    Qstart <- as.numeric(Qstart)
    Qend <- as.numeric(Qend)
    Sstart <- as.numeric(Sstart)
    Send <- as.numeric(Send)
    text_width <- as.numeric(text_width)
    Qstart_adjusted <- vector()
    Qtext_adjusted <- vector()
    Qend_adjusted <- vector()
    Sstart_adjusted <- vector()
    Send_adjusted <- vector()
    Stext_adjusted <- vector()
    middle_adjusted <- vector()
    line = ceiling(str_length(Qtext)/text_width)
    for(i in 0:(line - 1)){
        ##message("i is ", i)
        ##message("start is ", start)
        ##message("end is ", end)
        Qstart_adjusted[i+1] <- Qstart + i * text_width
        Qtext_adjusted[i+1] <- substr(Qtext, start, end)
        Qend_adjusted[i+1] <- Qstart_adjusted[i+1] + text_width -1
        ##Sstart_adjusted[i+1] <- Sstart+i*text_width
        Stext_adjusted[i+1] <- substr(Stext, start, end)
        ##Send_adjusted[i+1] <- Sstart_adjusted[i+1] + text_width -1
        middle_adjusted[i+1] <- substr(middle, start, end)
        if(Sstart <= Send){
            Sstart_adjusted[i+1] <- Sstart+i*text_width
            Send_adjusted[i+1] <- Sstart_adjusted[i+1] + text_width -1
        }else{
            Sstart_adjusted[i+1] <- Sstart - i*text_width
            Send_adjusted[i+1] <- Sstart_adjusted[i+1] - text_width +1
        }
        start = start + text_width
        end = start + text_width -1
    }
    Qend_adjusted[length(Qend_adjusted)] <- Qend
    Send_adjusted[length(Send_adjusted)] <- Send
    pmap(list(Qstart_adjusted,
              Qend_adjusted,
              Sstart_adjusted,
              Send_adjusted,
              Qtext_adjusted,
              Stext_adjusted,
              middle_adjusted), Alignment_singleLine) %>%
        paste(collapse="\n")

}

server <- function(input, output, session){
  
    ##custom_db <- c("LvTx")
    ##custom_db_path <- c("/LV_transcriptome/LvTx")
    ## Chunk for database selection
    observe({
        x <- input$dbtype
        if (x == "Ref Type: Genome"){
            y <- c("blastn", "tblastn")
        }else if(x== "Ref Type: Proteome"){
            y <- c("blastp")
        }else{}
    # Can also set the label and select items
        updateSelectInput(session, "program",
                          label = "Program",
                          choices = y)
      ##selected = tail(y, 1)
    })

    observe({
        x <- 11
        if (input$program != "blastn"){
            x <- 3
        }
        
        updateTextInput(session, "word_size",
                        label = "Word (W) length:",
                        value = x,
                        placeholder = x)
    })

    
    ## Indicate database selected
    db_selected <- reactive({
        if(input$Reference == "William 82 a2v1" && input$dbtype == "Ref Type: Genome"){
            paste0(db_path, db_set[1])
        }else if(input$Reference == "William 82 a2v1" && input$dbtype == "Ref Type: Proteome"){
            paste0(db_path, db_set[2])
        }else if(input$Reference == "W05" && input$dbtype == "Ref Type: Genome"){
            paste0(db_path, db_set[3])
        }else if(input$Reference == "William 82 a2v1" && input$dbtype == "Ref Type: Proteome"){
            paste0(db_path, db_set[4])
        }
    })
    
    blastresults <- eventReactive(input$blast, {
        ## Use shiny internal progress indicator instead
        withProgress(message = 'Running blast', value = 0, {
        ##gather input and set up temp file
            query <- input$query
            tmpInput <- tempfile(fileext = ".fa")
            tmpOutput <- tempfile(fileext = ".xml")
            #if else chooses the right database
            ##if (input$db == custom_db){
            ##  db <- db_selected()
            ##  remote <- c("")
            ##} else {
            ##  db <- c("nr")
            ##  #add remote option for nr since we don't have a local copy
            ##  remote <- c("-remote")
            ##}
            remote <- c("") # remote is disabled
            ## Disable multiple input sequence
            validate(
                need(str_count(query, ">") == 1, "Only one input is allowed.")
            )
            ## Check sequence
            query_stripped <- query %>% str_replace(., "^>.*\n", "")
            query_stripped <- str_replace_all(query_stripped, "\n", "")
            k <- str_split(query_stripped,"")[[1]]
            DNA_RNA_ALPHABET <- unique(c(Biostrings::DNA_ALPHABET, Biostrings::RNA_ALPHABET,
                                         tolower(Biostrings::DNA_ALPHABET), tolower(Biostrings::RNA_ALPHABET)))
            ##message("k is ", k)
            if (str_detect(input$dbtype,"Proteome")){
                validate(
                    need(sum(k %in% Biostrings::AA_ALPHABET) == length(k), "Please provide valid amino acid sequences.")
                )
            }else{
                validate(
                    need(sum(k %in% DNA_RNA_ALPHABET) == length(k), "Please provide valid DNA/RNA sequences.")
                )
            }
            ##this makes sure the fasta is formatted properly
            if (startsWith(query, ">")){
              writeLines(query, tmpInput)
            } else {
              writeLines(paste0(">Query\n",query), tmpInput)
            }
            
            #calls the blast
            system(paste0(input$program," -query ",tmpInput," -db ",db_selected()," -evalue ",input$eval," -word_size ",input$word_size," -num_alignments ",input$num_alignments," -outfmt 5 "," -out ",tmpOutput, remote))
            ##xmlParse(data)
            read_xml(tmpOutput, as = "text", encoding = "UTF-8")
        
        })
    }, ignoreNULL= T)

    ## First Extract program information
    ## This will be header of result table
    program_info <- reactive({
        if(is.null(blastresults())){}
        else {
            Extract_blastXMLheader(blastresults())
        }
    })


    output$program_version <- renderText({
        x <- program_info() %>% select(program_version) %>% pull(program_version)
        paste0("<b>Program:</b> ", x)
    })
    output$database_info <- renderText({
        x <- program_info() %>% select(db_name) %>% as.character()
        x <- x %>% str_split(simplify=TRUE, pattern="/") %>%
            unlist() %>% last()
        y <- program_info() %>% select(db_num) %>% as.character()
        z <- program_info() %>% select(db_len) %>% as.character()
        paste0("<b>Database:</b> ", x, " (", y, " sequences, ", z, " total letters)")
    })
    
    ## Now to parse the results...
    ## All relationship of Query and Subject will be parsed to nested tibble
    parsedresults <- reactive({
        if (is.null(blastresults())){}
        else{
            Parse_blastXML(blastresults())
        }
    })

    ## Makes the datatable
    output$blastResults <- renderDataTable({
        if (is.null(blastresults())){
        } else {
            blastXMLout <- parsedresults()
            blastXMLout %>%
                unnest(cols=c("Hsp_list")) %>%
                mutate(Identity=paste0(round(as.numeric(Hsp_identity)/as.numeric(Hsp_align_len), digits=2)*100,"%")) %>%
                mutate(Evalue=Hsp_evalue) %>%
                ##mutate(Evalue=round(as.numeric(Hsp_evalue), digits=5)) %>%
                select(query, subject, subject_len,
                       Identity,
                       Hsp_align_len,
                       Hsp_bit_score, Evalue) %>%
                dplyr::rename(Query=query,
                              Subject=subject,
                              Subject_len=subject_len,
                              Alignment_Len=Hsp_align_len,
                              Bitscore=Hsp_bit_score)
        }
    }, selection="single", rownames = FALSE)
  
    ##this chunk gets the alignemnt information from a clicked row
    output$selected_alignment <- renderText({
        if(is.null(input$blastResults_rows_selected)){}
        else{
            clicked = input$blastResults_rows_selected
            selected <- parsedresults() %>%
                unnest(cols=c("Hsp_list")) %>%
                dplyr::slice(clicked)
            query <- selected %>%
                select(query) %>% pull()
            score <- selected %>%
                select(Hsp_score) %>% pull()
            bitscore <- selected %>%
                select(Hsp_bit_score) %>% pull()
            evalue <- selected %>%
                select(Hsp_evalue) %>% pull()
            Hsp_identity <- selected %>%
                select(Hsp_identity) %>% pull()
            Hsp_align_len <- selected %>%
                select(Hsp_align_len) %>% pull()
            Hsp_gaps <- selected %>%
                select(Hsp_gaps) %>% pull() %>% as.numeric()
            Hit_frame <- selected %>%
                select(Hsp_hit_frame) %>% pull() %>% as.numeric()
            ##tibble(score, bitscore, evalue, Hsp_identity, Hsp_align_len, Hsp_gaps)
            identity <- as.numeric(Hsp_identity)/as.numeric(Hsp_align_len)
            identity <- round(identity, digits = 4)
            identity <- identity * 100
            identity <- paste0(identity, "%")

            gaps <- as.numeric(Hsp_gaps)/as.numeric(Hsp_align_len)
            gaps <- round(gaps, digits = 4)
            gaps <- gaps * 100
            gaps <- paste0(gaps, "%")
            if(Hit_frame > 0){
                strand = "Plus / Plus"
            }else{
                strand = "Plus / Minus"
            }
            paste0("Query: ", query, "<br/>", "Score = ", bitscore, " bits (", score, "),  Evalue = ",
                   evalue, "<br/>", "Identities = ",
                   Hsp_identity, "/", Hsp_align_len, " (", identity, "),  ",
                   " Gaps = ", Hsp_gaps, "/", Hsp_align_len, " (", gaps, ")",
                   "<br/>", "Strand = ", strand)
        }
    })
    
    output$alignment <- renderText({
        if(is.null(input$blastResults_rows_selected)){}
        else{
            clicked = input$blastResults_rows_selected
            selected <- parsedresults() %>%
                unnest(cols=c("Hsp_list")) %>%
                dplyr::slice(clicked)
            #tmp <- selected %>% select(Hsp_query_from)
            
            #print(selected["Hsp_query_from"])
            Qstart <- selected["Hsp_query_from"] %>% pull()
            Qend <- selected["Hsp_query_to"] %>% pull()
            Sstart <- selected["Hsp_hit_from"] %>% pull()
            Send <- selected["Hsp_hit_to"] %>% pull()
            Qtext <- selected["Hsp_qseq"] %>% pull() %>% as.character()
            Stext <- selected["Hsp_hseq"] %>% pull() %>% as.character()
            middle <- selected["Hsp_midline"] %>% pull()
            Adjust_alignment(Qstart, Qend, Sstart, Send, Qtext, Stext, middle, 80)
            #paste0("Query\t248\tATGGATCCGGAGGCTGCACGAACTGCTCGAGAATCTCTTGACCTGGCGTTCCATATGTCC\t307\n\t\t||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct\t1\tATGGATCCGGAGGCTGCACGAACTGCTCGAGAATCTCTTGACCTGGCGTTCCATATGTCC\t60\n")
        }
    })
##  #this chunk makes the alignments for clicked rows
##  output$alignment <- renderText({
##    if(is.null(input$blastResults_rows_selected)){}
##    else{
##      xmltop = xmlRoot(blastresults())
##      
##      clicked = input$blastResults_rows_selected
##        
##      #loop over the xml to get the alignments
##      align <- xpathApply(blastresults(), '//Iteration',function(row){
##        top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
##        mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
##        bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
##        rbind(top,mid,bottom)
##      })
##        
##      #split the alignments every 40 carachters to get a "wrapped look"
##      alignx <- do.call("cbind", align)
##      splits <- strsplit(gsub("(.{40})", "\\1,", alignx[1:3,clicked]),",")
##      
##      #paste them together with returns '\n' on the breaks
##      split_out <- lapply(1:length(splits[[1]]),function(i){
##        rbind(paste0("Q-",splits[[1]][i],"\n"),paste0("M-",splits[[2]][i],"\n"),paste0("H-",splits[[3]][i],"\n"))
##      })
##      unlist(split_out)
##    }
##  })
}
shinyApp(ui, server)
