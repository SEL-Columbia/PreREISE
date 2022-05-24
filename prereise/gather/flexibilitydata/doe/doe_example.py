import doe_data

download_path = "doe"
aggregated_path = "aggregated_doe_flexibility.csv"

doe_data.download_doe(download_path)
doe_data.aggregate_doe(download_path, aggregated_path)
