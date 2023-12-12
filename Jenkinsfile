pipeline {
    agent any
    stages {
        stage("install") {
            steps {
                sh 'mamba env update --file environment.yml --prefix ./env || mamba env create --force --file environment.yml --prefix ./env'
                sh 'echo conda environment: $CONDA_PREFIX'
            }
        }
        stage("unit tests") {
            steps {
                sh 'R devtools::test()'
            }
        }
    }
}