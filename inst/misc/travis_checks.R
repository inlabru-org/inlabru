system("echo '############################## RUNNING devtools::check() ##############################'")
check_result = devtools::check('/home/work', run_dont_test=FALSE, check_dir='/home')

system("echo '------------------------------ results ------------------------------------------------'")
print(str(check_result))

system("echo '------------------------------ checkdir contents --------------------------------------'")
list.files(check_result$checkdir)

system("echo '------------------------------ exit status --------------------------------------------'")
exit_status = length(check_result$errors) + length(check_result$warnings) + length(check_result$notes)
print(paste0("Exit status: ", exit_status))

# if (file.exists("/home/inlabru.Rcheck/tests/testthat.Rout.fail")){
#   system("echo '############################## TESTTHAT FAIL LOG ##############################'")
#   system("cat /home/inlabru.Rcheck/tests/testthat.Rout.fail")
# }

quit('no', status=exit_status)
