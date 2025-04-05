---
layout: default
title: 首页
---

# 欢迎来到我的计算化学学习博客

这里记录了我在冰表面吸附模拟和Python科学计算方面的学习经验和项目进展。

## 最新文章

<ul>
  {% for post in site.posts limit:5 %}
    <li>
      <h2><a href="{{ site.baseurl }}{{ post.url }}">{{ post.title }}</a></h2>
      <p>{{ page.date | date: "%Y-%m-%d" }}</p>
      <p>{{ post.excerpt }}</p>
    </li>
  {% endfor %}
</ul>
