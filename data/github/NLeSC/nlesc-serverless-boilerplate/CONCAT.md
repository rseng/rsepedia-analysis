# NLeSC serverless boilerplate

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3628237.svg)](https://doi.org/10.5281/zenodo.3628237)

## What is it?

The web application was bootstrapped with [Create React App](https://github.com/facebook/create-react-app).
The web application uses [AWS AppSync](https://aws.amazon.com/appsync/) and [Amplify](https://aws-amplify.github.io/docs/) to setup and run infrastructure.

It uses [AWS Batch](https://aws.amazon.com/batch/) to submit a Docker image as a job.

## Frontend

To have a interactive frontend the backend should first be setup.

### `yarn start`

Runs the app in the development mode.<br />
Open [http://localhost:3000](http://localhost:3000) to view it in the browser.

### `yarn test`

This will run tests of frontend code.

### `yarn build`

Builds the app for production to the `build` folder.<br />

## Backend

### Architecture

* Uses Graphql to communicate between web browser and server
* Uses AWS DynamoDB for persisting data
* Uses AWS S3 for storing input/output data files
* Uses AWS Cognito for authentication
* Uses AWS S3 for deployment
* Uses AWS S3 for hosting
* Uses AWS Batch job, written in Python, to run a computation in a Docker container.
* Uses AWS ECR to store Docker image for AWS Batch
* Uses AWS lambda to submit a AWS Batch job
* Uses AWS lambda to cancel a AWS Batch job
* Uses AWS lambda to listen for a AWS Batch job state changes
* Allows AWS Batch job to read/write to S3 and DynamoDB
* Allows AWS lambda functions to read/write to S3 and DynamoDB

[![Drawio diagram](docs/architecture.png)](docs/architecture.drawio)

### Installation

Requirements: 
* nodejs, tested with v12.13.1
* yarn, NodeJS package manager
* Docker, used for building batch job Docker image
* aws cli (pip install awscli)
* AWS account

Install amplify cli with

```sh
npm install -g @aws-amplify/cli@4.12.0
```

Amplify cli needs to be installed globally, to not pollute your env we suggest to use [nvm](https://github.com/nvm-sh/nvm) to isolate the node env.

Initialize amplify with

```sh
amplify configure
Follow these steps to set up access to your AWS account:

Sign in to your AWS administrator account:
https://console.aws.amazon.com/
Press Enter to continue

Specify the AWS Region
? region:  eu-central-1
Specify the username of the new IAM user:
? user name:  amplify-********
Complete the user creation using the AWS console
https://console.aws.amazon.com/iam/home?region=undefined#/users$new?step=final&accessKey&userNames=amplify-**********&permissionType=policies&policies=arn:aws:iam::aws:policy%2FAdministratorAccess
Press Enter to continue

Enter the access key of the newly created user:
? accessKeyId:   **********
? secretAccessKey:  **********
This would update/create the AWS Profile in your local machine
? Profile Name:  boiler

Successfully set up the new user.
```

To start when there are local amplify resources, but none in cloud with the amplify environment suffix.

```sh
amplify env add
# Asks for new environment name
```

To use current deployed amplify environment in the cloud

```sh
amplify env pull
```

To start when there are no local or cloud amplify resources (not the case for this boilerplate repo).

```sh
amplify init
```
The amplify resources that we made are logged in [docs/log.md](docs/log.md)


### To deploy services

Before amplify push you need to 
1. correct in subnet and security group in `amplify/backend/batch/task/parameters.json`.

```sh
# Deploy backend
amplify push
# Deploy frontend
amplify publish
```

The url where the application is running will be printed to screen.

After amplify you need to
1. Build & push the Docker image (see `amplify/backend/batch/task/` directory) which has the computation in it
2. Create users in [AWS Cognito console](https://eu-central-1.console.aws.amazon.com/cognito/home?region=eu-central-1).
### Add Docker image

The AWS Batch job definition has a Docker image. The Dockerfile is in `amplify/backend/batch/task/src` directory.

To build and push see View push commands popup on [AWS Container Repository](https://eu-central-1.console.aws.amazon.com/ecr/repositories/nlesc-hello-task-master/?region=eu-central-1)

For example
```
$(aws ecr get-login --no-include-email --region eu-central-1)
docker build -t nlesc-hello-task-master .
docker tag nlesc-hello-task-master:latest <account id>.dkr.ecr.eu-central-1.amazonaws.com/nlesc-hello-task-master:latest
docker push <account id>.dkr.ecr.eu-central-1.amazonaws.com/nlesc-hello-task-master:latest
```# Example graphql queries

## Create job

```graphql
mutation CreateJob(
  $input: CreateJobInput!
  $condition: ModelJobConditionInput
) {
  createJob(input: $input, condition: $condition) {
    id
  }
}
```

Vars
```json
{"input": {
    "status": {
        "state": "UNKNOWN",
        "submittedBy": "stefan",
        "submittedAt": "2020-01-24T10:03:54.337Z",
        "updatedAt": "2020-01-24T10:03:54.337Z"
    },
    "jobDescriptionId": "12be628b-ecb4-48e9-bd6a-9fbd01c64bf6",
    "owner": "stefan"
}

}
```# Log of Amplify commands ran to setup the backend

## Init amplify

```sh
amplify init
? Enter a name for the project nlesc
? Enter a name for the environment master
? Choose your default editor: Visual Studio Code
? Choose the type of app that you're building javascript
Please tell us about your project
? What javascript framework are you using react
? Source Directory Path:  src
? Distribution Directory Path: build
? Build Command:  yarn build
? Start Command: yarn start
Using default provider  awscloudformation

For more information on AWS Profiles, see:
https://docs.aws.amazon.com/cli/latest/userguide/cli-multiple-profiles.html

? Do you want to use an AWS profile? Yes
? Please choose the profile you want to use boiler

....

✔ Successfully created initial AWS cloud resources for deployments.
✔ Initialized provider successfully.
Initialized your environment successfully.

Your project has been successfully initialized and connected to the cloud!

Some next steps:
"amplify status" will show you what you've added already and if it's locally configured or deployed
"amplify <category> add" will allow you to add features like user login or a backend API
"amplify push" will build all your local backend resources and provision it in the cloud
“amplify console” to open the Amplify Console and view your project status
"amplify publish" will build all your local backend and frontend resources (if you have hosting category added) and provision it in the cloud

Pro tip:
Try "amplify add api" to create a backend API and then "amplify publish" to deploy everything
```

## Add API and database

```sh
amplify add api
? Please select from one of the below mentioned services: GraphQL
? Provide API name: nlesc
? Choose the default authorization type for the API Amazon Cognito User Pool
Using service: Cognito, provided by: awscloudformation
 
 The current configured provider is Amazon Cognito. 
 
 Do you want to use the default authentication and security configuration? Default configuration
 Warning: you will not be able to edit these selections. 
 How do you want users to be able to sign in? Username
 Do you want to configure advanced settings? No, I am done.
Successfully added auth resource
? Do you want to configure advanced settings for the GraphQL API No, I am done.
? Do you have an annotated GraphQL schema? No
? Do you want a guided schema creation? Yes
? What best describes your project: Single object with fields (e.g., “Todo” with ID, name, description)
? Do you want to edit the schema now? Yes
Please edit the file in your editor: .../nlesc-serverless-boilerplate/amplify/backend/api/nlesc/schema.graphql
? Press enter to continue 

The following types do not have '@auth' enabled. Consider using @auth with @model
	 - Todo
Learn more about @auth here: https://aws-amplify.github.io/docs/cli-toolchain/graphql#auth 


GraphQL schema compiled successfully.

Edit your schema at .../nlesc-serverless-boilerplate/amplify/backend/api/nlesc/schema.graphql or place .graphql files in a directory at .../nlesc-serverless-boilerplate/amplify/backend/api/nlesc/schema
Successfully added resource nlesc locally

Some next steps:
"amplify push" will build all your local backend resources and provision it in the cloud
"amplify publish" will build all your local backend and frontend resources (if you have hosting category added) and provision it in the cloud
```

So the Graphql schema is in [amplify/backend/api/nlesc/schema.graphql](amplify/backend/api/nlesc/schema.graphql).

## Push resource to AWS

```sh
amplify push
✔ Successfully pulled backend environment master from the cloud.

Current Environment: master

| Category | Resource name | Operation | Provider plugin   |
| -------- | ------------- | --------- | ----------------- |
| Auth     | nlesc8d53e119 | Create    | awscloudformation |
| Api      | nlesc         | Create    | awscloudformation |
? Are you sure you want to continue? Yes

The following types do not have '@auth' enabled. Consider using @auth with @model
	 - Todo
Learn more about @auth here: https://aws-amplify.github.io/docs/cli-toolchain/graphql#auth 


GraphQL schema compiled successfully.

Edit your schema at .../nlesc-serverless-boilerplate/amplify/backend/api/nlesc/schema.graphql or place .graphql files in a directory at .../nlesc-serverless-boilerplate/amplify/backend/api/nlesc/schema
? Do you want to generate code for your newly created GraphQL API Yes
? Choose the code generation language target typescript
? Enter the file name pattern of graphql queries, mutations and subscriptions src/graphql/**/*.ts
? Do you want to generate/update all possible GraphQL operations - queries, mutations and subscriptions Yes
? Enter maximum statement depth [increase from default if your schema is deeply nested] 5
? Enter the file name for the generated code src/API.ts

...

✔ Generated GraphQL operations successfully and saved at src/graphql
✔ Code generated successfully and saved in file src/API.ts
✔ All resources are updated in the cloud

GraphQL endpoint: https://******************.appsync-api.eu-central-1.amazonaws.com/graphql

```

### Create a cognito user & test Graphql endpoint

Without a web application with a registration form you can 
1. create an coginto user in [AWS console](https://eu-central-1.console.aws.amazon.com/cognito/home?region=eu-central-1)
2. Goto [AWS Appsync console](https://eu-central-1.console.aws.amazon.com/appsync/home?region=eu-central-1)
3. Goto Run a query
4. Login with the user you just created
5. Run following query

```graphql
{
  listTodos {
    items {
      name
    }
  }
}
```

Results in 

```json
{
  "data": {
    "listTodos": {
      "items": []
    }
  }
}
```

### Host the app

Create S3 bucket for hosting
```sh
amplify add hosting
? Select the environment setup: DEV (S3 only with HTTP)
? hosting bucket name nlesc-serverless-boilerplate
? index doc for the website index.html
? error doc for the website index.html

You can now publish your app using the following command:
Command: amplify publish
```

### Publish app

```sh
amplify publish
✔ Successfully pulled backend environment master from the cloud.

Current Environment: master

| Category | Resource name   | Operation | Provider plugin   |
| -------- | --------------- | --------- | ----------------- |
| Hosting  | S3AndCloudFront | Create    | awscloudformation |
| Auth     | nlesc8d53e119   | No Change | awscloudformation |
| Api      | nlesc           | No Change | awscloudformation |
? Are you sure you want to continue? Yes

...

Done in 37.05s.
frontend build command exited with code 0
✔ Uploaded files successfully.
Your app is published successfully.
http://nlesc-serverless-boilerplate-master.s3-website.eu-central-1.amazonaws.com
```

### Add markAsCompleted lambda function

Using Javascript

```sh
amplify function add
? Provide a friendly name for your resource to be used as a label for this category in the project: markAsCompleted
? Provide the AWS Lambda function name: nlesc-markAsCompleted
? Choose the function template that you want to use: Hello world function
? Do you want to access other resources created in this project from your Lambda function? Yes
? Select the category api
Api category has a resource called nlesc
? Select the operations you want to permit for nlesc read, update
```

### Add markAsCompleted mutation to Graphql

To `amplify/backend/api/nlesc/schema.graphql` added
```graphql
type Mutation {
  markAsCompleted(todoId: ID!): ID @function(name: "nlesc-markAsCompleted-${env}")
}
```

To `amplify/backend/function/markAsCompleted/markAsCompleted-cloudformation-template.json:Resources:AmplifyResourcesPolicy.Properties.PolicyDocument` add policy to allow dynomdb access.
```

Update Graphql api with
```sh
amplify push
```

### Upgrade node8 lambdas

Since Jan 2020 node8 lambdas have been deprecated (https://aws-amplify.github.io/docs/cli/lambda-node-version-update). 

```
npm install -g @aws-amplify/cli@4.12.0
amplify status
Amplify CLI uses Lambda backed custom resources with CloudFormation to manage part of your backend resources.
In response to the Lambda Runtime support deprecation schedule
https://docs.aws.amazon.com/lambda/latest/dg/runtime-support-policy.html
Nodejs runtime need to be updated from nodejs8.10  to nodejs10.x in the following template files:
/home/stefanv/git/marin/nlesc-serverless-boilerplate/amplify/backend/auth/nlesc8d53e119/nlesc8d53e119-cloudformation-template.yml

Please test the changes in a test environment before pushing these changes to production. There might be a need to update your Lambda function source code due to the NodeJS runtime update. Please take a look at https://aws-amplify.github.io/docs/cli/lambda-node-version-update for more information

? Confirm to update the NodeJS runtime version to 10.x Yes

NodeJS runtime version updated successfully to 10.x in all the CloudFormation templates.
Make sure the template changes are pushed to the cloud by "amplify push"

Current Environment: master

| Category | Resource name   | Operation | Provider plugin   |
| -------- | --------------- | --------- | ----------------- |
| Auth     | nlesc8d53e119   | Update    | awscloudformation |
| Api      | nlesc           | No Change | awscloudformation |
| Hosting  | S3AndCloudFront | No Change | awscloudformation |
| Function | markAsCompleted | No Change | awscloudformation |

GraphQL endpoint: https://rxcbtxh7gvc3naaqvgm5vcteka.appsync-api.eu-central-1.amazonaws.com/graphql
Hosting endpoint: http://nlesc-serverless-boilerplate-master.s3-website.eu-central-1.amazonaws.com


amplify push
```


### Add markAllAsCompleted lambda function

Will call `markAsCompleted` for all incomplete todos.

### Add markAllAsCompleted mutation to Graphql

### Add todoCounter lambda function

Using Python count the number of completed and incomplete todos.

### Add todoCounter query to Graphql

### Add Batch job definition

Create custom category to setup a AWS Batch compute environment, job queue and job definition.

```
cd amplify/backend
mkdir -p batch/task
```

Create cloud formation file `amplify/backend/batch/task/template.json` using https://aws.amazon.com/cloudformation/getting-started/ as a guide and https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-resource-batch-computeenvironment.html as a reference for the different resources. 

Batch job settings you might want to change in `amplify/backend/batch/task/template.json`:
* Resources:ComputeEnvironment:ComputeResources:InstanceTypes, defaults to optimal, but can be changed to for example only run `p2.family` instances to get only gpu vms.
* Desired number of vcpus (Resources:ComputeEnvironment:ComputeResources:DesiredvCpus), defaults to 0 which means no vm will be kept running when queue is empty.
* change Resources:BatchJobPolicy resource to give the Docker container access to other AWS resources
In `amplify/backend/batch/task/parameters.json`:
* taskImageName, Docker image name
* security group, goto [EC2 console](https://eu-central-1.console.aws.amazon.com/ec2/v2/home?region=eu-central-1#SecurityGroups:sort=groupId) to find your security group (adjust url to your region)
* subnets, goto [VPC console](https://eu-central-1.console.aws.amazon.com/vpc/home#subnets:sort=SubnetId) to find your subnets (adjust url to your region)
* taskMemory, memory in Mb for container
* taskVcpus, number of vCpus for container
* service role https://docs.aws.amazon.com/batch/latest/userguide/service_IAM_role.html
* instance role https://docs.aws.amazon.com/batch/latest/userguide/instance_IAM_role.html


Make cli aware of batch category
```
amplify env checkout master

Current Environment: master

| Category | Resource name   | Operation | Provider plugin   |
| -------- | --------------- | --------- | ----------------- |
| Batch    | task            | Create    | awscloudformation |
...
```

Create batch resources in cloud
```
amplify push
```

### Add Docker image

The AWS Batch job definition has a Docker image. The Dockerfile is in `amplify/backend/batch/task/src` directory.

To build and push see View push commands popup on [AWS Container Repository](https://eu-central-1.console.aws.amazon.com/ecr/repositories/nlesc-hello-task-master/?region=eu-central-1)

For example
```
$(aws ecr get-login --no-include-email --region eu-central-1)
docker build -t nlesc-hello-task-master .
docker tag nlesc-hello-task-master:latest <account id>.dkr.ecr.eu-central-1.amazonaws.com/nlesc-hello-task-master:latest
docker push <account id>.dkr.ecr.eu-central-1.amazonaws.com/nlesc-hello-task-master:latest
```

### Submit job from AWS console

1. Goto https://eu-central-1.console.aws.amazon.com/batch/home?region=eu-central-1#/jobs
2. Press `Submit job`
3. Fill job name
4. Select `nlesc-task-jobdefinition-master` as job definition
5. Select `nlesc-jobqueue-master` as job queue
6. Add Parameter called `jobdescriptionid`
7. Set vCpus==1 and memory==512
8. Press `Submit job`
9. In console follow progress
10. Once completed goto [Cloudwatch batch job logs](https://eu-central-1.console.aws.amazon.com/cloudwatch/home?region=eu-central-1#logStream:group=/aws/batch/job) to see the stdout/stderr

### Add job model to graphql

Push it
```
amplify push
```



### Adjust Docker image to talk to appsync/dynomdb

Job should fetch and update something from DynamoDB.

Add a job description in [AWS AppSync console](https://eu-central-1.console.aws.amazon.com/appsync/home?region=eu-central-1#/apis) with
```graphql
mutation CreateJobDescription(
  $input: CreateJobDescriptionInput!
) {
  createJobDescription(input: $input) {
    id
    payload {
      count
    }
  }
}
```
Variables:
```
{
  "input": {
    "payload": {
      "count": 42
    }
  }
}
```

Submit job from AWS Batch console using job description id returned by mutation.

### Add job submit function lambda

```
amplify function add
? Provide a friendly name for your resource to be used as a label for this category in the project: jobsubmit
? Provide the AWS Lambda function name: nlesc-jobsubmit
? Choose the function template that you want to use: Hello world function
? Do you want to access other resources created in this project from your Lambda function? Yes
? Select the category api
Api category has a resource called nlesc
? Select the operations you want to permit for nlesc create, read, update
```

After function body and permission have been set push again
```
amplify push
```

### Add job submit mutation to graphql

Add mutation to graphql.
Test with AWS AppSync console

Graphql:
```graphql
mutation SubmitJob($jobdescriptionid: ID!) {
  submitJob(jobdescriptionid: $jobdescriptionid)
}
```
Vars:
```
{
  "jobdescriptionid": "b78be18c-24e0-4fc2-8837-9ef4495d7cff"
}
```
Should return a job id and create a job in the db.

### Add job cancel function lambda

```
amplify function add
Using service: Lambda, provided by: awscloudformation
? Provide a friendly name for your resource to be used as a label for this category in the project: jobcancel
? Provide the AWS Lambda function name: nlesc-jobcancel
? Choose the function template that you want to use: Hello world function
? Do you want to access other resources created in this project from your Lambda function? Yes
? Select the category api
Api category has a resource called nlesc
? Select the operations you want to permit for nlesc read
```


After function body and permission have been set push again
```
amplify push
```

### Addd progress reporting to Docker

Use boto3 to update DynamoDB job.status.

### Add job cancel mutation to graphql

Add mutation to graphql.
Test with AWS AppSync console with the following query
```
 mutation CancelJob {
   cancelJob(jobid: "4b57218f-6364-49f3-a742-340a82593bd4")
 }
```
Check in AWS Batch console if job has been cancelled

### Add job listen function lambda

Added lambda function and attached it to the AWS batch job state change event.

After function body and permission have been set push again
```
amplify push
```

Submit a new job and check in DynamoDB if job status has been updated (aka not SUBMITTED).


### Add authorization to Graphql models and functions/docker

Add auth and owner fields to graphql.

Check React app still works 
Check job submissions from console still works.

### Add job functionality to React app

* [x] CRUD for job descriptions
* [x] List of jobs with status and result
* [x] Button to cancel active job
* [x] submit job using existing job description
