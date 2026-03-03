#' Fill effects based on individuals state and time spent in state
#'
#' @param M_t   health states occupied by individuals at time t
#' @param cl    cycle length
#' @param k     simulation number (or iteration)
#'
#'
#' @note Additional information on variables used in function:
#'   - m_PRE: matrix tracking time spent in PRE state
#'   - m_REL1: matrix tracking time spent in REL1 state
#'   - m_REL2: matrix tracking time spent in REL2 state
#'   - m_BMT1: matrix tracking time spent in BMT1 state
#'   - m_BMT1_REL: matrix tracking time spent in BMT1_REL state
#'   - j: current time/cycle index
#'   - u_PRE: utility values for PRE state (vector, indexed by k)
#'   - u_PRE_LT: utility values for PRE state long-term (vector, indexed by k)
#'   - u_REL: utility values for REL states (vector, indexed by k)
#'   - u_REL_LT: utility values for REL states long-term (vector, indexed by k)
#'   - u_BMT: utility values for BMT1 state (vector, indexed by k)
#'   - u_BMT_LT: utility values for BMT1 state long-term (vector, indexed by k)
#'   - u_DEAD: utility value for DEAD
#' @note Defining these utility values is done in '01_input_params.R'
#'
#' @return vector of QALYs for individuals during cycle t
#'
Effs <- function(M_t, cl = 30 / 365, k = k) {
  u_t <- rep(0, length(M_t)) # By default the utility for everyone is zero

  u_t[M_t == "PRE"  & m_PRE[,j+1]  <  365*2] <- u_PRE[k]          # Post-diagnosis (no relapse hx) short-term (if <2 years in state) [Furlong et al., Table 3 HUI3]
  u_t[M_t == "PRE"  & m_PRE[,j+1]  >= 365*2] <- u_PRE_LT[k]       # Post-diagnosis (no relapse hx) long-term  (if >2 years in state) [Furlong et al., Table 3 HUI3]
  u_t[M_t == "REL1" & m_REL1[,j+1] <  365*5] <- u_REL[k]          # Post-relapse short-term (if <5 years in state) [Kelly et al. Table 1, EQ-5D]
  u_t[M_t == "REL1" & m_REL1[,j+1] >= 365*5] <- u_REL_LT[k]       # Post-relapse long-term (if >5 years in state) [Lawitschka et al, Table 3, PedsQL]
  u_t[M_t == "REL2" & m_REL2[,j+1] <  365*5] <- u_REL[k]          # Post-2nd relapse short-term (if <5 years in state) [Kelly et al. Table 1, EQ-5D]
  u_t[M_t == "REL2" & m_REL2[,j+1] >= 365*5] <- u_REL_LT [k]      # Post-2nd relapse short-term (if <5 years in state) [Kelly et al. Table 1, EQ-5D]
  u_t[M_t == "BMT1" & m_BMT1[,j+1] <  365*5] <- u_BMT[k]          # Post-BMT short-term (if <5 years in state) [Kurosawa et al., Table 5, EuroQOL5D]
  u_t[M_t == "BMT1" & m_BMT1[,j+1] >= 365*5] <- u_BMT_LT[k]       # Post-BMT long-term (if >5 years in state)  [Lawitschka et al., Table 3, PedsQL]
  u_t[M_t == "BMT1_REL" & m_BMT1_REL[,j+1] <  365*5] <- u_REL[k]  # Relapse post-BMT short-term (if <5 years in state) [Kelly et al. Table 1, EQ-5D]
  u_t[M_t == "BMT1_REL" & m_BMT1_REL[,j+1] >= 365*5] <- u_REL[k]  # Relapse post-BMT long-term (if >5 years in state)  [Lawitschka et al., Table 3, PedsQL]
  u_t[M_t == "DEAD"] <- u_DEAD

  QALYs <- u_t * cl # calculate the QALYs during cycle t
  return(QALYs)
}
