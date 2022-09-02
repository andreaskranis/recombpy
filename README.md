
Install the project as ```pip install -e <path-to-rootDir>```

In __init__.py I use this approach: https://towardsdatascience.com/whats-init-for-me-d70a312da583
So, you need to ensure that the src.<packageName>.__init__.py has the right imports

For **testing**, check [this guide](https://realpython.com/python-testing/). I am using nose2

For **documentation** I use the [google style](https://stackoverflow.com/questions/3898572/what-is-the-standard-python-docstring-format)

=======
# recombpy
 
from macbook 13