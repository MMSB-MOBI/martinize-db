# Martinize Database Server

## Create JWT keys

```bash
mkdir .keys
cd .keys
# do not enter passphrase
ssh-keygen -t rsa -b 4096 -m PEM -f key_new.pem
openssl rsa -in key_new.pem -pubout -outform PEM -out key_new
```
