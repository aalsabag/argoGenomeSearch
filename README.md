
### Installing Argo
```bash
kubectl create namespace argo
kubectl apply -n argo -f https://raw.githubusercontent.com/argoproj/argo/stable/manifests/install.yaml
```
I recommend also getting the CLI which you can find for your respective operating system [here](https://github.com/argoproj/argo/releases).
For windows you will have to download the exe and rename it to `argo.exe` and add it to your path.