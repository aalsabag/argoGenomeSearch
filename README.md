
## Prerequisites
1. [Argo](#argo)
2. [Minio or a similar artifact storage tool](#minio)

## Getting Started

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
$ helm install argo-artifacts stable/minio --set service.type=LoadBalancer --set fullnameOverride=argo-artifa
```
And then it should be running on port 9000 by default!

#### If you would instead like a dedicated running minio instance to the drive of your choosing:
Download from [here](https://min.io/download#/)
And then add it to your path and run it like so:
```
minio server F:/Data #or whatever drive you want to use
```
Ideally this would be remote.
