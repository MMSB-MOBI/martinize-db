martinize_db_client_dir="../martinize-db-client/build"

rm -rf ./static/www/*
cp -R $martinize_db_client_dir/* ./static/www

mkdir deploy

cp -R .keys build static templates package.json settings.json deploy/

zip -r deploy.zip deploy

rm -rf deploy
