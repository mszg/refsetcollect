from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='home'),
    path('run/', views.run_sampling, name='run_sampling'),
    path('progress/<str:job_id>/', views.progress, name='progress'),
    path('result/<str:job_id>/', views.result, name='result'),
    path('download-genome/', views.download_genome_ncbi, name='download_genome'),
    path('download-selected/<str:job_id>/', views.download_selected_genomes, name='download_selected_genomes'),
]
