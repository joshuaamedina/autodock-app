all: deploy job submit

deploy:
	tapis apps deploy -W ./

job:
	tapis jobs init --no-notify --no-archive --output vina_job.json python_vina-0.0.1

submit:
	tapis jobs submit -F vina_job.json

history:
	tapis jobs history ${JOB}
