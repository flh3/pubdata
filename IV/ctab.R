ctab <- function(data, row_var, col_var, percent_type = "row",
                 show_na = TRUE) {
   
   require(janitor)
   require(dplyr)
   
   data %>%
      tabyl({{ row_var }}, {{ col_var }}, show_na = show_na) %>%
      adorn_totals(c("row", "col")) %>%
      adorn_percentages(percent_type) %>%
      adorn_pct_formatting(digits = 1) %>%
      adorn_ns(position = "front") %>%
      adorn_title("combined")
   
}