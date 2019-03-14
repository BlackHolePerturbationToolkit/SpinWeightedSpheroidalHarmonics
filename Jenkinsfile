pipeline {
  agent {
    docker {
      image 'wolfram-docker-11.3.0 '
    }

  }
  stages {
    stage('Run tests') {
      agent {
        docker {
          image 'wolfram-docker-11.3.0'
        }

      }
      steps {
        sh 'Tests/AllTests.wls'
        junit 'TestReport.xml'
      }
    }
  }
}