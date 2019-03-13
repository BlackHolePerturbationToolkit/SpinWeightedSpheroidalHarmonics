pipeline {
  agent any
  stages {
    stage('Run tests') {
      agent any
      steps {
        sh 'Tests/AllTests.wls'
        junit 'TestReport.xml'
      }
    }
  }
}