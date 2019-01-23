system("echo '############################## RUNNING CRAN CHECK ##############################'")
cran_checks = do.call(c, as.list(devtools::check('/home/work', run_dont_test=FALSE, check_dir='/home')))

if (file.exists("/home/inlabru.Rcheck/tests/testthat.Rout.fail")){
  system("echo '############################## TESTTHAT FAIL LOG ##############################'")
  system("cat /home/inlabru.Rcheck/tests/testthat.Rout.fail")
}

quit('no', status=length(cran_checks))
