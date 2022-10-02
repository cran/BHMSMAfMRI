.onAttach <- function(libname, pkgname) {
  message <- c("\n Welcome! Thanks for trying BHMSMAfMRI.",
               "\n \n Website: https://nilotpalsanyal.github.io/BHMSMAfMRI/",
               "\n Bug report: https://github.com/nilotpalsanyal/BHMSMAfMRI/issues")
  packageStartupMessage(message)
}