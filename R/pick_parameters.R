
#-------------------------------------------------------------------------------

pick_parameters <- function(
    bp_rule = "runif(1, min = 0.5, max = 3)",
    dp_rule,
    rq_rule = NULL,
    bq_rule = "0",
    dq_rule = "0",
    t_start_dist,
    t_run_rule,
    nr_of_families = NULL,
    nr_burst_divs = 3,
    correction_t_burst = NULL,
    response_nr = 1,
    prev_parameters = NULL,
    quality_dist = NULL,
    ASD = FALSE) {

  parameters <- data.frame(cell_type = factor(),
                           div_counter = integer(),
                           t_start = numeric(),
                           quality = numeric(),
                           fam_nr = integer(),
                           fam_nr_2 = integer(),
                           fam_nr_3 = integer(),
                           nr_burst_divs = integer(),
                           t_burst = numeric(),
                           t_run = numeric(),
                           t_correct = numeric(),
                           bp = numeric(),
                           dp = numeric(),
                           bq = numeric(),
                           dq = numeric())

  if (is.null(prev_parameters) == FALSE) {
    prev_q <- prev_parameters[prev_parameters$cell_type == "Q", ]
    nr_of_families <- nrow(prev_q)
    prev_divs <- prev_q$div_counter
  } else {
    prev_divs <- rep(0, nr_of_families)
  }

  if (!is.null(quality_dist)) {
    quality <- round(eval(parse(text = paste(quality_dist))), digits = 2)
    # beta distribution: rbeta(nr_of_families, shape1 = 2, shape2 = 5)
    # uniform distribution: runif(nr_of_families, 0, 1)
  } else {
    quality <- rep(NA, nr_of_families)
  }

  nr_burst_divs <- eval(parse(text = nr_burst_divs))
  div_counter <- nr_burst_divs + prev_divs
  t_burst <- rep(0.2, nr_of_families)
  t_start_expr <- substitute(eval(parse(text = t_start_dist)))
  t_start_val <- eval(t_start_expr)
  t_start <- round(t_start_val, digits = 2)
  max_run_time <- 8

  for (i in 1:nr_of_families) {
    # Make sure the families have the correct family number
    if (response_nr == 1) {
      fam_nr <- i
      fam_nr_2 <- 0
      fam_nr_3 <- 0
    }
    if (response_nr == 2) {
      fam_nr <- prev_q$fam_nr[i]
      fam_nr_2 <- i
      fam_nr_3 <- 0
    }
    if (response_nr == 3) {
      fam_nr <- prev_q$fam_nr[i]
      fam_nr_2 <- prev_q$fam_nr_2[i]
      fam_nr_3 <- i
    }

    frac_rq <- eval(parse(text = rq_rule))
    q_cells <- ifelse(ASD == FALSE,
                      min(
                          max(ceiling(frac_rq * 2^nr_burst_divs[i]), 0),
                          2^nr_burst_divs[i]),
                      1)
    p_cells <- (2^nr_burst_divs[i]) - q_cells

    if (p_cells != 0) {
      for (cell in 1:p_cells) {
        # First define the expression and then evaluate the expression
        # If this is not included, the expression is only evaluated once,
        # leading to the same t_run for all (sub)families
        t_run_expr <- substitute(eval(parse(text = t_run_rule)))
        t_run_val <- eval(t_run_expr)
        t_run <- round(t_run_val, digits = 2)

        # this is to make sure all proliferation stops at t=7
        if (t_start[i] + t_run >= max_run_time) {
          t_correct <- round((t_start[i] + t_run - max_run_time), digits = 2)
        } else {
          t_correct <- 0
        }
        if (t_start[i] < 0) {
          t_start[i] <- 0.01
        }
        if (t_start[i] > 6.5) {
          t_start[i] <- 6.5
        }
        if (t_start[i] < 2) {
          t_start[i] <- 2.5
        }

        bp <- round(eval(parse(text = bp_rule)), digits = 2)
        dp <- round(eval(parse(text = dp_rule)), digits = 2)

        parameters %<>% add_row(
          cell_type = "P",
          div_counter = div_counter[i],
          t_start = t_start[i],
          quality = quality[i],
          fam_nr = fam_nr,
          fam_nr_2 = fam_nr_2,
          fam_nr_3 = fam_nr_3,
          nr_burst_divs = nr_burst_divs[i],
          t_burst = t_burst[i],
          t_run = t_run,
          t_correct = t_correct,
          bp = bp,
          dp = dp,
          bq = 0,
          dq = 0
        )
      }
    }

    if (q_cells != 0) {
      for (cell in 1:q_cells) {
        t_run <- 0
        # this is to make sure all proliferation stops at t=7
        if (t_start[i] + t_run >= max_run_time) {
          t_correct <- round((-max_run_time + t_start[i] + t_run), digits = 2)
        } else {
          t_correct <- 0
        }
        if (t_start[i] < 0) {
          t_start[i] <- 0.01
        }
        if (t_start[i] > 6.5) {
          t_start[i] <- 6.5
        }

        bq <- round(eval(parse(text = bq_rule)), digits = 2)
        dq <- round(eval(parse(text = dq_rule)), digits = 2)

        parameters %<>% add_row(
          cell_type = "Q",
          div_counter = div_counter[i],
          t_start = t_start[i],
          quality = quality[i],
          fam_nr = fam_nr,
          fam_nr_2 = fam_nr_2,
          fam_nr_3 = fam_nr_3,
          nr_burst_divs = nr_burst_divs[i],
          t_burst = t_burst[i],
          t_run = 0,
          t_correct = t_correct,
          bp = 0,
          dp = 0,
          bq = bq,
          dq = dq
        )
      }
    }
  }
  return(parameters)
}
