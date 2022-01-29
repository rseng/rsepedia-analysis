<div align="center">
  <br>
  <h1>QPixel</h1>
  <strong>Q&A by the community, for the community</strong>
</div>
<br>
<p align="center">
  <a href="https://circleci.com/gh/codidact/qpixel">
    <img src="https://circleci.com/gh/codidact/qpixel.svg?style=svg" alt="CircleCI Build Status">
  </a>
  <a href="https://coveralls.io/github/codidact/qpixel">
    <img src="https://coveralls.io/repos/github/codidact/qpixel/badge.svg" alt="Coverage Status">
  </a>
  <a href="https://zenodo.org/badge/latestdoi/237078806">
    <img src="https://zenodo.org/badge/237078806.svg" alt="DOI">
  </a>
</p>

Rails-based version of our core software, powering [codidact.com](https://codidact.com). Currently under active development towards MVP.

## Table of Contents
- [Community](#community)
- [Contributing](#contributing)
- [Installation](#installation)

## Community
To discuss features and voice ideas, please ask a new question on [Codidact Meta](https://meta.codidact.com). For technical discussions about the QPixel software itself, please ask on [Codidact Collab](https://collab.codidact.org) instead.

## Contributing
Contributions are welcome - please read the [CONTRIBUTING](https://github.com/codidact/qpixel/blob/develop/CONTRIBUTING.md)
document before you start and look at the [GitHub issues](https://github.com/codidact/qpixel/issues) for things to do.

## Installation
These instructions are assuming you already have a Unix environment available with Ruby and Bundler installed.
WSL should work as well, but (core) Windows is unlikely to.

If you don't already have Ruby installed, use [RVM](https://rvm.io/) or
[rbenv](https://github.com/rbenv/rbenv#installation) to install it before following these instructions.

### Install prerequisites

For Debian-Based Linux:

```
sudo apt update
sudo apt install gcc
sudo apt install make
sudo apt install libmysqlclient-dev
sudo apt install autoconf bison build-essential libssl-dev libyaml-dev libreadline-dev zlib1g-dev libncurses5-dev libffi-dev libgdbm-dev
sudo apt install mysql-server
```

For Arch-Based Linux:

```
sudo pacman -Syyu
sudo pacman -Sy gcc
sudo pacman -Sy make
sudo pacman -Sy ruby autoconf bison base-devel unixodbc
sudo pacman -Sy openssl
sudo pacman -S mariadb mysqld nodejs
```

For Mac:

```
xcode-select --install
brew install mysql bison openssl mysql-client
bundle config --global build.mysql2 --with-opt-dir="$(brew --prefix openssl)"
```

### 1. Install JS runtime
If you already have Node.JS installed, you can skip this step. If not,
[download and install it](https://nodejs.org/en/download/).

### 2. Install Redis
If you haven't already got it, [download and install Redis](https://redis.io/download).

### 3. Install Imagemagick

If you haven't already installed Imagemagick, you'll need to [install it for
your system](https://imagemagick.org/script/download.php).

### 4. Download QPixel
Clone the repository and `cd` into the directory:

    git clone https://github.com/codidact/qpixel
    cd qpixel

### 5. Configure database connection
If you weren't asked to set the root MySQL user password during `mysql-server` installation, the installation is
likely to be using Unix authentication instead. You'll need to sign into the MySQL server with `sudo mysql -u root`
and create a new database user for QPixel:

```sql
CREATE USER qpixel@localhost IDENTIFIED BY 'choose_a_password_here';
GRANT ALL ON qpixel_dev.* TO qpixel@localhost;
GRANT ALL ON qpixel_test.* TO qpixel@localhost;
GRANT ALL ON qpixel.* TO qpixel@localhost;
```

Copy `config/database.sample.yml` to `config/database.yml` and fill in the correct host, username, and password for
your environment. If you've followed these instructions (i.e. you have installed MySQL locally), the correct host
is `localhost` or `127.0.0.1`.

You'll also need to fill in details for the Redis connection. If you've followed these instructions, the sample file
should already contain the correct values for you, but if you've customised your setup you'll need to correct them.

### 6. Set up QPixel
Install gems:

    bundle install

Set up the database:

    rails db:create
    rails db:schema:load
    rails r db/scripts/create_tags_path_view.rb
    rails db:migrate

 You'll need to create a Community record and purge the Rails cache before you can seed the database. In a Rails
 console (`rails c`), run:

```ruby
Community.create(name: 'Dev Community', host: 'localhost:3000')
Rails.cache.clear
```

### 7. Seed the database:

    $ rails db:seed
    Category: Created 2, skipped 0
    [...]

Run the server!

    rails s
    
### 8. Create admin account    
    
You can create the first user account in the application, which should be running at `http://localhost:3000/`. To upgrade
the user account
to an admin account, run `rails c` for a console, followed by:

```ruby
User.last.update(confirmed_at: DateTime.now, is_global_admin: true)
```    

If you create more accounts, you can visit `http://localhost:3000/letter_opener` to see the emails that would otherwise be sent by QPixel.

### 9. Configure Categories

Before you try to create a post we need to configure categories! 
Go to `http://localhost:3000/categories/`

![img/categories.png](img/categories.png)

 Click "edit" for each category and scroll down to see the "Tag Set" field. This
 will be empty on first setup.

![img/tagset.png](img/tagset.png)

You will need to select a tag set for each category! For example, the Meta category can be
associated with the "Meta" tag set, and the Q&A category can be associated with "Main"

![img/tagset-selected.png](img/tagset-selected.png)

Make sure to click save for each one.<br> 
<em>Note:</em> You may need to run `rails db:seed` again.

### 10. Create a Post

You should then be able to create a post! There are character requirements for the
body and title, and you are required at least one tag.

![img/create-post.png](img/create-post.png)

And then click to "Save Post in Q&A"

![img/post.png](img/post.png)


### Install with Docker

See the README.md in the [docker](docker) folder for complete instructions.

## License
[AGPL licensed](https://github.com/codidact/qpixel/blob/master/LICENSE).

<br>

[⬆ Back to Top](#table-of-contents)
# Security Policy

## Supported Versions

While this project is under development (i.e. until release 1.0.0), **only the current version is supported**. If you
believe you have found a security vulnerability, ensure you are running the latest commit from GitHub; if you're not,
update to it &mdash; if you are, please report the vulnerability to us.

## Reporting a Vulnerability

Please [report a security issue](https://codidact.atlassian.net/servicedesk/customer/portal/1/group/1/create/17) via our
support portal. You will receive an automated confirmation that we have received your report, and we'll aim to make an
initial assessment and get back to you within 24 hours. This timeframe is not a guarantee, and some cases may take longer.

In all cases, please *do not* disclose any details of any security issue, potential or confirmed, until we have confirmed
to you that it has been either fixed or confirmed not a vulnerability.

## Bug Bounty & Acknowledgements
There is no bug bounty program available for this software at present. If you wish us to make an acknowledgement of your
contribution, please let us know under what name you'd like it to be added, and to what URL you would like us to link.
# Contributing
Contributing to QPixel follows broadly the same process as any other Codidact project.

## What needs doing?
 - Most bugs and change requests are here on GitHub. Have a look at the open issues to find something that needs doing.
 - There are a few other items in the [TODO list in the wiki](https://github.com/codidact/qpixel/wiki/TODO-list).
   
Once you've picked what you're going to work on, please **leave a comment** on the issue to indicate you're planning to work on
it; this helps us reduce wasted effort. If there's not already an issue for the feature you want to work on, please create one.
If you need time to work on an issue, that's absolutely fine, but please **keep us updated** with comments on the issue - if we
don't hear from you for a few weeks, we may assume you've given up working on that issue and give it to someone else.

## What's the workflow?
 * First, **you need an issue to work under**. Either pick an existing issue or create a new one, and leave a comment on it
   to indicate that you're working on it.
 * Second, you can make your changes. If you have write access to the repository, create a topic branch (please use the format
   `art/40/add-bells-and-whistles`, i.e. `username/issue-number/brief-description`) and make your changes there; if not, fork
   the repository and work in your fork.
 * Once you've made your changes, submit a pull request targeting the `develop` branch.

Keep in mind that **status checks are required to pass** and **at least one approving review** is required from the team before
any pull request can be merged. If status checks don't pass, we won't be able to merge - there are no exceptions, so please fix
the failures and commit again. You can always mark your pull request as a draft while you're still trying to make it work.

## What standards are there?
We have code style and standards documents for each applicable language. Please make sure you follow these if possible; if
there's a good reason why not, please document it in your code, add a linter exception, and let us know why in your pull
request. Here they are:

 * [Code standards: CSS](https://github.com/codidact/core/wiki/Code-standards:-CSS)
 * [Code standards: CSS naming](https://github.com/codidact/core/wiki/Code-standards:-CSS-naming)
 * [Code standards: HTML](https://github.com/codidact/core/wiki/Code-standards:-HTML)
 * [Code standards: JS](https://github.com/codidact/core/wiki/Code-standards:-JS)
 * There is a .rubocop.yml file provided in the project and rubocop is included in the bundle; please run `bundle exec rubocop` for 
   Ruby style checking.
 
When writing CSS, keep in mind that our design framework, [Co-Design](https://design.codidact.org/) is available in QPixel, and
should be used where possible. Avoid writing custom CSS if you can; favour using components and atomic classes from Co-Design.

We also have some [guidelines for commit messages](https://github.com/codidact/core/wiki/Committing-guidelines). Again, please
follow these where possible, as they help us to keep a cohesive commit history and see how the project has developed.
# Docker Installation

A [docker-compose.yml](../docker-compose.yml) file is provided for deployment with Docker compose, if you choose.

## 1. Build Containers

You should first build the images, before making any changes to config files. We do this so the container
is not built with secrets.

```bash
docker-compose build
```

If you need to just rebuild one container, you can do that too.

```bash
docker-compose build uwsgi
docker-compose build db
docker-compose build redis
```

## 2. Setup and Secrets

The `docker-compose.yml` file uses a `.env` file in the same directory to load dynamic values used when the docker containers are initialized.

This is useful for setting up custom values locally. Additionally, your secrets (the mysql database credentials and admin user name) are inserted into the running container through the `docker/env` file.

Both the `.env` file and the `docker/env` file are gitignored, so you can change values to suit, these filed need to be copied to the correct locations with some default values. You can do this in one step by executing a bash script.

```bash
# ensure script is executable, from the project root:
chmod +x docker/local-setup.sh
docker/local-setup.sh
```

Editing the `./.env` file will modify the corresponding variables used in the docker-compose.yml file but **NOT** the environment variables in the container. Editing the `./docker/env` file will change environment variables only in the running container.

## 3. Database File
Ensure `config/database.yml` has the username and password as defined in [docker/env](dummy.env) file. The `config/database.yml` should already be gitignored.

The `COMMUNITY_NAME` value defined in the `.env` file defines the initial community name on your local DB.

the `COMMUNITY_ADMIN_USERNAME`, `COMMUNITY_ADMIN_PASSWORD` and `COMMUNITY_ADMIN_EMAIL` values in the `docker/env` file define the first user you can log in as - however you will need to follow the instructions below to ensure you can log in as that user.

## 4. Start Containers

Then start your containers! 

```bash
docker-compose up # append -d if you want to detach the processes, although it can be useful to see output into the terminal
Creating qpixel_redis_1 ... done
Creating qpixel_db_1    ... done
Creating qpixel_uwsgi_1 ... done
```

The uwsgi container has a sleep command for 15 seconds to give the database a chance to start,
so don't expect to see output right away. After about 20 seconds, check to make sure the server is running (and verify port 3000, note that you can change this mapping in the `.env` file)

```
uwsgi_1  | => Booting Puma
uwsgi_1  | => Rails 5.2.4.3 application starting in development 
uwsgi_1  | => Run `rails server -h` for more startup options
uwsgi_1  | Puma starting in single mode...
uwsgi_1  | * Version 3.12.6 (ruby 2.6.5-p114), codename: Llamas in Pajamas
uwsgi_1  | * Min threads: 0, max threads: 16
uwsgi_1  | * Environment: development
uwsgi_1  | * Listening on tcp://localhost:3000
uwsgi_1  | Use Ctrl-C to stop
```

You should then be able to open your browser to [http://localhost:3000](http://localhost:3000)
and see the interface. 

![img/interface.png](../img/interface.png)

Before you login, since we don't have email configured, you'll need to set a manual
`confirmed_at` variable for your newly created user. You can do this easily with a single
command to the container:

```bash
$ docker exec qpixel_uwsgi_1 rails runner "User.second.update(confirmed_at: DateTime.now)"
Running via Spring preloader in process 111
```

The first user is the system user, so the second user is the admin created during the
start of the container. And you can of course do this same command for any future users that you don't want to require email confirmation for. You can then click "Sign in" to login with what you defined for `$COMMUNITY_ADMIN_EMAIL` and `$COMMUNITY_ADMIN_PASSWORD`. Importantly, your password must be 6 characters or greater, otherwise the user won't be created. 

## 5. Login

Once you are logged in, you should see your icon in the top right:

![img/logged-in.png](../img/logged-in.png)

## 6. Configure Categories

Before you try to create a post we need to configure categories! 
Go to `http://localhost:3000/categories/`

![img/categories.png](../img/categories.png)

 Click "edit" for each category and scroll down to see the "Tag Set" field. This
 will be empty on first setup.

![img/tagset.png](../img/tagset.png)

You will need to select a tag set for each category! For example, the Meta category can be
associated with the "Meta" tag set, and the Q&A category can be associated with "Main"

![img/tagset-selected.png](../img/tagset-selected.png)

Make sure to click save for each one.

## 7. Create a Post

You should then be able to create a post! There are character requirements for the
body and title, and you are required at least one tag.

![img/create-post.png](../img/create-post.png)

And then click to "Save Post in Q&A"

![img/post.png](../img/post.png)

That's it!

### 8. Stop Containers

When you are finished, don't forget to clean up.

```bash
docker-compose stop
docker-compose rm
```

### 9. Next steps

The current goal of this container is to provide a development environment for
working on QPixel. This deployment has not been tested with email notifications
enabled / set up, nor to deploy in production mode. If you require these 
modifications, please [open an issue](https://github.com/codidact/qpixel/issues).
We've noticed that you've had some recent interactions on this site that have been less than polite.

In keeping with our [Code of Conduct](/policy/code-of-conduct), we expect everyone using this site to remain respectful, presume good intent, and use common sense when communicating with other members of the site.

We understand that a single person is rarely solely to blame when discussions get heated, and that you may feel that you are not at fault. However, please remember that someone else being out of line does not excuse similar behavior on your part. If someone is being disrespectful or rude, please raise a flag on the comment or post and then disengage from the interaction. A moderator will handle your flag and review the situation.

While we appreciate your continued contributions to $SiteName, we do ask that you stay mindful of your fellow users and abide by the Code of Conduct.We've noticed that in addition to your contributions to $SiteName, you've also been using this account to promote a product, service, or similar.

This is a gentle reminder that using an account for promotional purposes is against the rules of this site.

Please take a moment to review the guidelines for promotional content here on $SiteName, which can be found in the [help center](/policy/spam).

While we appreciate your continued contributions to $SiteName, we do ask that you stay within the boundaries of the rules about promotional content.We've noticed that you've written a number of posts about topics that are beyond the scope of this site.

You can find out about what's on-topic and what's off-topic on $SiteName in the [help center](/help/faq).

This is just a gentle reminder that we expect posts on this site to stay focused on the topic on-hand.
We have a [Network of communities](https://codidact.com/) that you are free to use; you may find one of our other communities more suitable to some of your posts. If you would like to talk about the possibility of creating a site for a subject not currently covered, please write a post on [meta](https://meta.codidact.com/categories/10) with your request.
Additionally, we have a $ChatLink available for more free-form discussion.

While we appreciate your continued contributions within the scope of this site, we do ask that you make sure that the topic of your posts remain in scope.