system("echo '############################## RUNNING CRAN CHECK ##############################'")
cran_checks = do.call(c, as.list(devtools::check('/home/work', run_dont_test=FALSE, check_dir='/home')))

system("echo '############################## TESTTHAT LOG ##############################'")
system("cat /home/inlabru.Rcheck/tests/testthat.Rout.fail")

quit('no', status=length(cran_checks))
