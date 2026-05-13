import click
import pandas as pd

@click.command()
@click.option('-i', "--input", required=True, help='File for list of coverages.')
@click.option('-e', "--empty", required=True, help="File to be written is no genome > min coverage")
@click.option('-c', "--cutoff", default=0.1, help='Minimum % genome coverage.', show_default=True)
def check_empty(input, empty, cutoff):
    cov = pd.read_table(input, usecols=['gotu', 'coverage_ratio'])
    print(cov['coverage_ratio'])
    print(cutoff)
    print((cov['coverage_ratio'] > cutoff))
    print((cov['coverage_ratio'] > cutoff).sum())
    if (cov['coverage_ratio'] > cutoff).sum():
        with open(empty, 'w') as o:
            pass


if __name__ == "__main__":
    check_empty()
