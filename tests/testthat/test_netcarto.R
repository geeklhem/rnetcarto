library(rnetcarto)

context("Test the result of the R interface 'netcarto'.")
test_that("netcarto works with a simple input matrix", {
    input = matrix(0,6,6)
    input[1,2] = 1
    input[2,3] = 1
    input[3,1] = 1
    input[3,4] = 1
    input[4,5] = 1
    input[5,6] = 1
    input[6,4] = 1
    rownames(input) = c("A","B","C","D","E","F")
    colnames(input) = rownames(input)
    ans = netcarto(input)
    expect_equal(ans[[2]],0.3571429,tolerance=0.01)
    expect_equal(length(ans),2)
    expect_equal(nrow(ans[[1]]),6)
    expect_equal(ncol(ans[[1]]),4)
})

test_that("netcarto works with a simple input list", {
    nd1 = c("A","B","C","D","E","F","C")
    nd2 = c("B","C","A","E","F","D","D")
    ans = netcarto(list(nd1,nd2))
    expect_equal(ans[[2]],0.3571429,tolerance=0.01)
    expect_equal(length(ans),2)
    expect_equal(nrow(ans[[1]]),6)
    expect_equal(ncol(ans[[1]]),4)
})
