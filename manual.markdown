---
layout: page
title: Manual
permalink: /manual/
order: 4
---

{% for file in site.static_files %}
    {% if file.path contains '/doc/functions/' and file.extname == '.html' %}
- [{{ file.name | remove: ".html" }}]({{ site.baseurl }}{{ file.path }})
    {% endif %}
{% endfor %}