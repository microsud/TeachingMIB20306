plot_volcano_custom <-
  function (dep,
            contrast,
            label_size = 3,
            add_names = TRUE,
            adjusted = FALSE)
  {
    if (is.integer(label_size))
      label_size <- as.numeric(label_size)
    assertthat::assert_that(
      inherits(dep, "SummarizedExperiment"),
      is.character(contrast),
      length(contrast) == 1,
      is.numeric(label_size),
      length(label_size) == 1,
      is.logical(add_names),
      length(add_names) ==
        1,
      is.logical(adjusted),
      length(adjusted) == 1
    )
    row_data <- rowData(dep)
    if (any(!c("name", "ID") %in% colnames(row_data))) {
      stop(
        paste0(
          "'name' and/or 'ID' columns are not present in '",
          deparse(substitute(dep)),
          "'.\nRun make_unique() to obtain required columns."
        ),
        call. = FALSE
      )
    }
    if (length(grep("_p.adj|_diff", colnames(row_data))) < 1) {
      stop(
        paste0(
          "'[contrast]_diff' and '[contrast]_p.adj' columns are not present in '",
          deparse(substitute(dep)),
          "'.\nRun test_diff() to obtain the required columns."
        ),
        call. = FALSE
      )
    }
    if (length(grep("_significant", colnames(row_data))) < 1) {
      stop(
        paste0(
          "'[contrast]_significant' columns are not present in '",
          deparse(substitute(dep)),
          "'.\nRun add_rejections() to obtain the required columns."
        ),
        call. = FALSE
      )
    }
    if (length(grep(paste(contrast, "_diff", sep = ""), colnames(row_data))) ==
        0) {
      valid_cntrsts <-
        row_data %>% data.frame() %>% select(ends_with("_diff")) %>%
        colnames(.) %>% gsub("_diff", "", .)
      valid_cntrsts_msg <- paste0("Valid contrasts are: '",
                                  paste0(valid_cntrsts, collapse = "', '"),
                                  "'")
      stop(
        "Not a valid contrast, please run `plot_volcano()` with a valid contrast as argument\n",
        valid_cntrsts_msg,
        call. = FALSE
      )
    }
    diff <-
      grep(paste(contrast, "_diff", sep = ""), colnames(row_data))
    if (adjusted) {
      p_values <- grep(paste(contrast, "_p.adj", sep = ""),
                       colnames(row_data))
    }
    else {
      p_values <- grep(paste(contrast, "_p.val", sep = ""),
                       colnames(row_data))
    }
    signif <- grep(paste(contrast, "_significant", sep = ""),
                   colnames(row_data))
    df <- data.frame(
      x = row_data[, diff],
      y = -log10(row_data[,
                          p_values]),
      z = row_data[, signif],
      name = row_data$name
    )
    signif_prots <- df %>% filter(z)
    name1 <- gsub("_vs_.*", "", contrast)
    name2 <- gsub(".*_vs_", "", contrast)
    p <- ggplot(df, aes(x, y)) + geom_vline(xintercept = 0) +
      geom_point(col = "grey") + geom_point(data = signif_prots,
                                            col = "steelblue") + geom_text(data = data.frame(),
                                                                       aes(
                                                                         x = c(Inf,-Inf),
                                                                         y = c(-Inf,-Inf),
                                                                         hjust = c(1, 0),
                                                                         vjust = c(-1,-1),
                                                                         label = c(name1, name2),
                                                                         size = 5,
                                                                         fontface = "bold"
                                                                       )) +
      labs(x = "Fold change (log2)") + theme_DEP1() + theme(legend.position = "none")
    if (add_names) {
      p <- p + ggrepel::geom_text_repel(
        data = signif_prots,
        aes(label = name),
        size = label_size,
        box.padding = unit(0.1,
                           "lines"),
        point.padding = unit(0.1, "lines"),
        segment.size = 0.2
      )
    }
    if (adjusted) {
      p <- p + labs(y = "Adjusted p-value (-log10)")
    }
    else {
      p <- p + labs(y = "P-value (-log10)")
    }
    p
  }
