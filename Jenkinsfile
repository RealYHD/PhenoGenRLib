pipeline {
    agent any
    stages {
        stage("install") {
            steps {
                sh 'mamba env update --file environment.yml --prefix ./env || mamba env create --force --file environment.yml --prefix ./env'
                sh 'Rscript scripts/install.R'
                sh 'echo conda environment: $CONDA_PREFIX'
            }
        }
        stage("check") {
            steps {
                sh 'Rscript scripts/check.R'
            }
        }
        stage("unit tests") {
            steps {
                sh 'Rscript scripts/test.R'
                xunit checksName: '', tools: [JUnit(excludesPattern: '', pattern: 'tests/testthat/test_results.xml', stopProcessingIfError: true)]
            }
        }
        stage("build") {
            steps {
                sh 'Rscript scripts/build.R'
            }
        }
    }
}