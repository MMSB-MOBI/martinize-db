# Configure Docker for computers with proxies

## Mac & Windows

Server don't work on Windows. Please don't use Windows.

On Mac OS, install Docker Desktop. Configure proxies in Docker Desktop GUI, in Dashboard > Settings > Resources > Proxies section.


## Linux 

It assumes you use Debian or one of its derivates, and your system service organizer is `systemd`.

### Install docker
```bash
sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common

# Replace debian by ubuntu if you use ubuntu
curl -fsSL https://download.docker.com/linux/debian/gpg | sudo apt-key add -

### DEBIAN
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/debian \
   $(lsb_release -cs) \
   stable"

### UBUNTU
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"


sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io

# Docker should be installed.

sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker 

mkdir ~/.docker

# Enable docker service on startup
sudo systemctl enable docker
sudo systemctl restart docker
```

### Install docker-compose

Ensure $HTTP_PROXY, $http_proxy, $HTTPS_PROXY and $https_proxy variables are set.

```bash
sudo curl -L "https://github.com/docker/compose/releases/download/1.25.4/docker-compose-Linux-x86_64" -o /usr/local/bin/docker-compose
sudo chmod +x /usr/local/bin/docker-compose
sudo ln -s /usr/local/bin/docker-compose /usr/bin/docker-compose
```

### Configure proxy settings

Create `~/.docker/config.json` file with the following content:

```json
{
 "proxies":
 {
   "default":
   {
     "httpProxy": "<insert here proxy>",
     "httpsProxy": "<insert here proxy>",
     "ftpProxy": "<insert here proxy>"
   }
 }
}
```

Configure the systemd proxy settings.

```bash
sudo mkdir -p /etc/systemd/system/docker.service.d
printf "[Service] \n\
Environment=\"HTTP_PROXY=<INSERT HERE THE PROXY>\"" > /etc/systemd/system/docker.service.d/http-proxy.conf

printf "[Service] \n\
Environment=\"HTTPS_PROXY=<INSERT HERE THE PROXY>\"" > /etc/systemd/system/docker.service.d/https-proxy.conf

sudo systemctl daemon-reload

sudo systemctl restart docker

# Check if the proxy is set
docker info | grep Proxy
```

## Check if everything is okay

Run the following command:

```bash
docker run hello-world
```
