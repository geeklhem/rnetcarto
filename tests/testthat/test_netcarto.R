library(rnetcarto)

context("Test simple non bipatite networks.")

test_that("Two triplets", {
    input = matrix(0,6,6)
    input[1,2] = 1
    input[2,3] = 1
    input[3,1] = 1
    input[4,5] = 1
    input[5,6] = 1
    input[6,4] = 1
    rownames(input) = c("A","B","C","D","E","F")
    colnames(input) = rownames(input)
    ans = netcarto(input)
    #print(ans)
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[2], label="Module 1 (A/B)")
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[3], label="Module 1 (A/C)")
    expect_equal(ans[[1]]$'module'[2],ans[[1]]$'module'[3], label="Module 1 (B/C)")

    expect_equal(ans[[1]]$'module'[4],ans[[1]]$'module'[5], label="Module 2 (A/B)")
    expect_equal(ans[[1]]$'module'[4],ans[[1]]$'module'[6], label="Module 2 (A/C)")
    expect_equal(ans[[1]]$'module'[5],ans[[1]]$'module'[6], label="Module 2 (B/C)")

    for (i in 1:6) {
      expect_equal(ans[[1]]$'participation'[i],0, label="participation")
      expect_equal(ans[[1]]$'connectivity'[i],0, label="participation")
    }
    expect_equal(ans[[2]],0,tolerance=0.01,label="Modularity")
})
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
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[2], label="Module 1 (A/B)")
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[3], label="Module 1 (A/C)")
    expect_equal(ans[[1]]$'module'[2],ans[[1]]$'module'[3], label="Module 1 (B/C)")

    expect_equal(ans[[1]]$'module'[4],ans[[1]]$'module'[5], label="Module 2 (A/B)")
    expect_equal(ans[[1]]$'module'[4],ans[[1]]$'module'[6], label="Module 2 (A/C)")
    expect_equal(ans[[1]]$'module'[5],ans[[1]]$'module'[6], label="Module 2 (B/C)")

    expect_equal(ans[[2]],0.3571429,tolerance=0.01,label="Modularity")
    expect_equal(length(ans),2)
    expect_equal(nrow(ans[[1]]),6)
    expect_equal(ncol(ans[[1]]),4)
})

test_that("netcarto works with a simple input list", {
    nd1 = c("A","B","C","D","E","F","C")
    nd2 = c("B","C","A","E","F","D","D")
    ans = netcarto(list(nd1,nd2))
    expect_equal(ans[[2]],0.3571429,tolerance=0.01,label="Modularity")
    expect_equal(length(ans),2)
    expect_equal(nrow(ans[[1]]),6)
    expect_equal(ncol(ans[[1]]),4)
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[2], label="Module 1 (A/B)")
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[3], label="Module 1 (A/C)")
    expect_equal(ans[[1]]$'module'[2],ans[[1]]$'module'[3], label="Module 1 (B/C)")

    expect_equal(ans[[1]]$'module'[4],ans[[1]]$'module'[5], label="Module 2 (A/B)")
    expect_equal(ans[[1]]$'module'[4],ans[[1]]$'module'[6], label="Module 2 (A/C)")
    expect_equal(ans[[1]]$'module'[5],ans[[1]]$'module'[6], label="Module 2 (B/C)")
})

context("Test simple bipatite networks.")
test_that("netcarto works with a simple bipartite input matrix", {
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
    ans = netcarto(input,seed=1,bipartite=TRUE)

    expect_equal(length(ans),2)
    expect_equal(nrow(ans[[1]]),7)
    expect_equal(ncol(ans[[1]]),4)
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[2], label="Module 1 (A/B)")
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[3], label="Module 1 (A/C)")
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[4], label="Module 1 (A/D)")
    expect_equal(ans[[1]]$'module'[5],ans[[1]]$'module'[6], label="Module 2 (E/F)")
    expect_equal(ans[[1]]$'module'[5],ans[[1]]$'module'[7], label="Module 2 (E/G)")
    expect_equal(ans[[1]]$'connectivity',
                 tolerance=0.01,
                 c(0.57735027, -1.73205081,  0.57735027,  0.57735027,0,0,0),
                 label="connectivity")
    expect_equal(ans[[1]]$'participation',
                 c(0,0,0,.46875,0.4444,0.4444,0.4444),
                 tolerance=0.01,
                 label="Participation Coefficient")
})

test_that("netcarto works with a simple bipartite input list", {
    nd1 = c("A","B","C","D","A","C","D","D","E","F","G")
    nd2 = c("a","a","a","a","b","b","b","c","c","c","c")

    ans = netcarto(list(nd1,nd2),seed=1,bipartite=TRUE)

    expect_equal(length(ans),2)
    expect_equal(nrow(ans[[1]]),7)
    expect_equal(ncol(ans[[1]]),4)
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[2], label="Module 1 (A/B)")
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[3], label="Module 1 (A/C)")
    expect_equal(ans[[1]]$'module'[1],ans[[1]]$'module'[4], label="Module 1 (A/D)")
    expect_equal(ans[[1]]$'module'[5],ans[[1]]$'module'[6], label="Module 2 (E/F)")
    expect_equal(ans[[1]]$'module'[5],ans[[1]]$'module'[7], label="Module 2 (E/G)")
    expect_equal(ans[[1]]$'connectivity',
                 tolerance=0.01,
                 c(0.57735027, -1.73205081,  0.57735027,  0.57735027,0,0,0),
                 label="connectivity")
    expect_equal(ans[[1]]$'participation',
                 c(0,0,0,.46875,0.4444,0.4444,0.4444),
                 tolerance=0.01,
                 label="Participation Coefficient")

})
