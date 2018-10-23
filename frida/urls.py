from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('calculate/', views.calculate, name='calculate'),
    path('doc/params1/',views.docs_p1,name='docs-p1'),
    path('doc/params2/',views.docs_p2,name='docs-p2'),
    path('doc/params3/',views.docs_p3,name='docs-p3')
] 
