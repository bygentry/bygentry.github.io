---
layout: page
title: Analysis Portfolio
---

## Data Analysis in R
    bygentry
    last update: {% last_modified_at %}

### AN588
#### [Homework 1](https://bygentry.github.io/AN588/HTMLs/homework1.html)
#### [Homework 2](https://bygentry.github.io/AN588/HTMLs/homework2.html)
#### [Homework 3](https://bygentry.github.io/AN588/HTMLs/homework3.html)
#### Homework 4
#### [Homework 5](https://bygentry.github.io/AN588/HTMLs/homework5.html)
#### Discriminant Analysis Module
#### Analysis Recreation Project


[Bootstrap4Pages Documentation](https://nicolas-van.github.io/bootstrap-4-github-pages/).



## License

[See the license file.](./LICENSE.md)

## Customization Guide

### Modify the configuration

You should at least edit the `_config.yml` file to edit your website's metadata, like the title, description and repository URL.

### Change your theme

This website uses the [Minty](https://bootswatch.com/minty/) Bootstrap theme by default. And you don't want to use the same theme everyone else uses do you?

You can of course modify anything in the `_includes`, `_layouts` and `_sass` folders to customize both the HTML or CSS of your website, possibly referring to the [Bootstrap documentation](https://getbootstrap.com/) or the [Jekyll documentation](https://jekyllrb.com/) when needed. This is a normal part of web development and it is outside the scope of this guide.

But if you don't know where to start I can recommend you to import a theme from [Bootswatch](https://bootswatch.com/).

* Go on [Bootswatch](https://bootswatch.com/) and choose a theme that you like.
* Using the top bar, download its `_variables.scss` and `_bootswatch.scss` files.
* Copy the content of `_variables.scss` in `_sass/_variables.scss`.
* Copy the content of `_bootswatch.scss` in `_sass/_bootstrap_customization.scss`.

That's it, you now have a totally different appearance for you website.

### Modify the content

You probably don't want the present guide to be the front page of your website, so you should edit the `index.md` file. You probably also want to edit or delete the `CONTRIBUTING.md`, `README.md` and `LICENSE.md` files.

Aside from that you can of course create new pages and posts like with any Jekyll website by refering to the [Jekyll documentation](https://jekyllrb.com/).

### Run Jekyll on your computer to speed up testing

Editing your website's content or theme directly on Github is completely possible but, due to the time Github Pages takes to update your website, it will probably be much more effective to work using a local Jekyll installation.

To do so:

* Install the [requirements for Jekyll](https://jekyllrb.com/docs/installation/).
* Type `bundle install` at the root of your project to install the necessary Ruby dependencies.
* Type `bundle exec jekyll serve` to launch the test Jekyll web server that will re-compile your work if you edit it.
* You can then open `http://localhost:4000` in your web browser to see your work-in-progress website.

Please note that, to ensure maximum compatibility with Github Pages, the `Gemfile` of this project references the `github-pages` gem, not Jekyll directly. This implies some differences in behavior compared to the official documentation of Jekyll.

## Known issues

* Bootstrap 5 should normally be post-processed using [Autoprefixer](https://github.com/postcss/autoprefixer). Even if it is possible to use autoprefixer with Jekyll, it is not possible with a classic Github Pages installation without adding some kind of pre-processing before publication. Since this project mostly aims compatibility with Github Pages I prefer to keep it that way. The consequences of this choice is that some Bootstrap features could not work as expected on older browsers.

## Other Github Pages related projects

I'm a fan of Github Pages for the possibilities it offers to anyone to publish a website for free. I have multiple projects that could be of interest if that's your case too:

* [Easy Markdown to Github Pages](https://nicolas-van.github.io/easy-markdown-to-github-pages/)
* [Parcel Github Pages Boilerplate](https://github.com/nicolas-van/parcel-github-pages-boilerplate)

