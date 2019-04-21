import os
import pytest

from .test_yaehmop import *

def test():
    pytest.main([os.path.dirname(__file__)])
