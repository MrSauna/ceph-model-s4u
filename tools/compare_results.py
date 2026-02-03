client_metrics = snakemake.input.client_metrics
mon_metrics = snakemake.input.mon_metrics
net_metrics = snakemake.input.net_metrics

for i in zip(client_metrics, mon_metrics, net_metrics):
    print("moi")
    print(i)