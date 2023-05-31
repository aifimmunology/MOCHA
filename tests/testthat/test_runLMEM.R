if (requireNamespace("lmerTest", quietly=TRUE)){
  
  test_that("runLMEM errors when failing to converge on a 1 sample test dataset", {
    cellPopulations <- c("C2", "C5")
  
    capture.output(
      ExperimentObj <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResults,
        cellPopulations = cellPopulations,
        threshold = 0
      ),
      type = "message"
    )
    
    modelFormula <- exp ~ Sample+PassQC
    
    expect_error(
      modelList <- runLMEM(
        ExperimentObj,
        modelFormula,
        initialSampling = 5,
        verbose = FALSE,
        numCores = 1
      ),
      "For the initial sampling, every test model failed to converge"
    )
  })
  
  test_that("runLMEM errors when failing to converge on a 3-sample dataset", {
    cellPopulations <- c("C2", "C3")
    capture.output(
      ExperimentObj <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResultsMultisample,
        cellPopulations = cellPopulations,
        threshold = 0
      )
    )
  
    modelFormula <- exp ~ Sample+PassQC
    
    expect_error(
      modelList <- runLMEM(
        ExperimentObj,
        modelFormula,
        initialSampling = 5,
        verbose = FALSE,
        numCores = 1
      ),
      "For the initial sampling, every test model failed to converge"
    )
  })
  
  test_that("pilotLMEM errors with formula on a 1-sample dataset", {
    cellPopulations <- c("C2", "C5")
  
    capture.output(
      ExperimentObj <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResults,
        cellPopulations = cellPopulations,
        threshold = 0
      ),
      type = "message"
    )
    
    modelFormula <- exp ~ Sample+PassQC
    
    expect_error(modelList <- pilotLMEM(
      ExperimentObj,
      cellPopulation = "C2",
      pilotIndices = c(1:5),
      modelFormula = exp ~ Sample+PassQC
    ), "No random effects terms specified in formula")
  
  })
  
  test_that("pilotLMEM errors with formula on a 3-sample dataset", {
    cellPopulations <- c("C2", "C3")
    capture.output(
      ExperimentObj <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResultsMultisample,
        cellPopulations = cellPopulations,
        threshold = 0
      )
    )
    
    expect_error(modelList <- pilotLMEM(
      ExperimentObj,
      cellPopulation = "C2",
      pilotIndices = c(1:5),
      modelFormula = exp ~ Sample+PassQC
    ), "No random effects terms specified in formula")
  
  })
  
  test_that("runLMEM errors with invalid formula", {
    cellPopulations <- c("C2", "C3")
    capture.output(
      ExperimentObj <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResultsMultisample,
        cellPopulations = cellPopulations,
        threshold = 0
      )
    )
  
    expect_error(modelList <- runLMEM(
      ExperimentObj,
      modelFormula = Sample ~ PassQC,
      initialSampling = 5,
      numCores = 1
    ), "modelFormula must start with 'exp' as the response.")
    
    expect_error(modelList <- runLMEM(
      ExperimentObj,
      modelFormula = exp ~ PassQC+Nonsense,
      initialSampling = 5,
      numCores = 1
    ), "Model factors are not found in the 'colData' of the ExperimentObj")
    
    expect_error(modelList <- runLMEM(
      ExperimentObj,
      modelFormula = "exp ~ PassQC+Sample",
      initialSampling = 5,
      numCores = 1
    ), "modelFormula is not a formula")
  })
  
  test_that("pilotLMEM errors with invalid formula", {
    cellPopulations <- c("C2", "C3")
    capture.output(
      ExperimentObj <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResultsMultisample,
        cellPopulations = cellPopulations,
        threshold = 0
      )
    )
  
    expect_error(modelList <- pilotLMEM(
      ExperimentObj,
      modelFormula = Sample ~ PassQC,
      cellPopulation = "C2",
    ), "modelFormula must start with 'exp' as the response.")
    
    expect_error(modelList <- pilotLMEM(
      ExperimentObj,
      modelFormula = exp ~ PassQC+Nonsense,
      cellPopulation = "C2",
    ), "Model factors are not found in the 'colData' of the ExperimentObj")
    
    expect_error(modelList <- pilotLMEM(
      ExperimentObj,
      modelFormula = "exp ~ PassQC+Sample",
      cellPopulation = "C2",
    ), "modelFormula is not a formula")
  })
}