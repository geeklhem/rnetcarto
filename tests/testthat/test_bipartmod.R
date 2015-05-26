library(rnetcarto)

context("Test the result of the R interface 'bipartmod'.")
test_that("bipartmod works with a simple input matrix", {
    input = matrix(0,7,3)
    input[1,1] = 1
    input[1,2] = 1
    input[2,1] = 1
    input[3,1] = 1
    input[3,2] = 1
    input[4,1] = 1
    input[4,2] = 1
    input[4,3] = 1
    input[5,3] = 1
    input[6,3] = 1
    input[7,3] = 1
    rownames(input) = c("A","B","C","D","E","F","G")
    colnames(input) = c("a",'b','c')
    ans = bipartmod(input,seed=1)

    expect_equal(length(ans),2)
    expect_equal(nrow(ans[[1]]),7)
    expect_equal(ncol(ans[[1]]),4)
    expect_equal(ans[[1]]$'module',c(2,2,2,2,1,1,1), label="Modules")
    expect_equal(ans[[1]]$'z-score',
                 tolerance=0.01,
                 c(0.57735027, -1.73205081,  0.57735027,  0.57735027,0,0,0),
                 label="Z-score")
    expect_equal(ans[[1]]$'participation',
                 c(0,0,0,.46875,0.4444,0.4444,0.4444),
                 tolerance=0.01,
                 label="Participation Coefficient")
})

test_that("bipartmod works with a simple input list", {
    nd1 = c("A","B","C","D","A","C","D","D","E","F","G")
    nd2 = c("a","a","a","a","b","b","b","c","c","c","c")

    ans = bipartmod(list(nd1,nd2),seed=1)

    expect_equal(length(ans),2)
    expect_equal(nrow(ans[[1]]),7)
    expect_equal(ncol(ans[[1]]),4)
    expect_equal(ans[[1]]$'module',c(2,2,2,2,1,1,1), label="Modules")
    expect_equal(ans[[1]]$'z-score',
                 tolerance=0.01,
                 c(0.57735027, -1.73205081,  0.57735027,  0.57735027,0,0,0),
                 label="Z-score")
    expect_equal(ans[[1]]$'participation',
                 c(0,0,0,.46875,0.4444,0.4444,0.4444),
                 tolerance=0.01,
                 label="Participation Coefficient")

})
