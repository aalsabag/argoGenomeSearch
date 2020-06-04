# Argo Genome Search
The Argo genome search project is an extension of the k8s genome search project. It's aim is to use argo to execute a pipeline that eventually runs parallel jobs to search for a particular sequence of amino acids. This is in contrast to the creation of plain k8s jobs executing in parallel
## Prerequisites
1. [Argo](#argo)
2. [Minio or a similar artifact storage tool](#minio)

![execution](./execution.gif)

## Getting Started
We begin by creating a container that will split a given faa file in proper FASTA format. The logic for splitting the files can be found in [split.py](./split.py). It takes in a file name as a parameter as well as a number defining how many files it should be split into. To create the docker container we execute:
```
docker build . -f Dockerfile.split -t filesplit:latest
```
Next we create a container that will search a file for us. The logic for searching a file can be found in [search.py](./search.py). It takes in a file as a parameter as well as a search string defining a sequence of amino acids.
```
docker build . -f Dockerfile.search -t filesearch:latest
```
Now we can add a FASTA formatted FAA file in our current working directory. I have one called influenze_xsmall.faa

Now we can submit our workflow through the argo cli. 
```
argo submit workflow.yaml --parameter sequence=NDVTSLISTTYPYTGPPPMSHGSSTKYT --parameter file=influenza_xsmall.yaml --parameter split-file-count=10
```
Here we specify the sequence we are searching for, the name of the file in our local dir, and the number of parallel search jobs we want to have running. To use default parameters, we can just run `argo submit workflow.yaml`

[workflow.yaml](./workflow.yaml) uses volume mounts as a method for compiling all of the results in a single file. This is until this [issue](https://github.com/argoproj/argo/issues/934) is resolved. You can see the reasons why this might be a bad idea depending on what type of infrastructure you are using. Otherwise, a basic method with uncompiled results (results in different files) can be found in [workflow-no-compile.yaml](./workflow-no-compile.yaml).
### [Installing Argo](#argo)
```bash
kubectl create namespace argo
kubectl apply -n argo -f https://raw.githubusercontent.com/argoproj/argo/stable/manifests/install.yaml
```
I recommend also getting the CLI which you can find for your respective operating system [here](https://github.com/argoproj/argo/releases).
For windows you will have to download the exe and rename it to `argo.exe` and add it to your path.

### [Installing Minio (Our artifact storage location)](#minio)
A helm chart has been created to configure and use minio so all you need to do is:
1. Make sure you have helm 3 installed (`choco install kubernetes-helm`)
2. Execute the following.
```
$ helm repo add stable https://kubernetes-charts.storage.googleapis.com/ # official Helm stable charts
$ helm repo update
$ helm install argo-artifacts stable/minio --set service.type=LoadBalancer --set fullnameOverride=argo-artifacts
```
And then it should be running on port 9000 by default!
3. Create a bucket called `my-bucket` . The default creds are
```
AccessKey: AKIAIOSFODNN7EXAMPLE
SecretKey: wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY
```
4. Configure the workflow controller to use minio
```
kubectl edit cm -n argo workflow-controller-configmap
```
```
data:
  artifactRepository: |
    s3:
      bucket: my-bucket
      endpoint: argo-artifacts.default:9000
      insecure: true
      # accessKeySecret and secretKeySecret are secret selectors.
      # It references the k8s secret named 'argo-artifacts'
      # which was created during the minio helm install. The keys,
      # 'accesskey' and 'secretkey', inside that secret are where the
      # actual minio credentials are stored.
      accessKeySecret:
        name: argo-artifacts
        key: accesskey
      secretKeySecret:
        name: argo-artifacts
        key: secretkey
```
4. Access the argo UI
```
kubectl -n argo port-forward deployment/argo-server 2746:2746
```

#### If you would instead like a dedicated running minio instance to the drive of your choosing:
Download from [here](https://min.io/download#/)
And then add it to your path and run it like so:
```
minio server F:/Data #or whatever drive you want to use
```
Ideally this would be a remote server.
