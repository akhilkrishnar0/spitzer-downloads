from spitzer_pipeline_class import SpitzerPipeline

if __name__ == '__main__':
    pipe = SpitzerPipeline(galaxy_centers={'NGC1566': (65.0042, -54.9378)})
    pipe.run(['NGC1566'], bands=['3.6','4.5','5.8','8.0','24'],
             radius_arcmin=6, skip_download=True, do_plot=True)
