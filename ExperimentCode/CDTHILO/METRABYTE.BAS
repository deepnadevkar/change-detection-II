Attribute VB_Name = "Module3"
Option Explicit

'MetraByte
'Global Const METRABYTE_BASE_ADDRESS = &H308
Global Const METRABYTE_BASE_ADDRESS = &H300
Global Const METRABYTE_PORT_B = METRABYTE_BASE_ADDRESS + 1
Global Const METRABYTE_CONTROL = METRABYTE_BASE_ADDRESS + 3
Global Const APPARATUS_BOGUS_LIGHT = &H40
Global Const APPARATUS_JUICE = &H20
Global Const APPARATUS_PELLETS = &H10
Global Const APPARATUS_HOUSE_LIGHT = &H8
Global Const APPARATUS_SET_B_OUTPUT = &H80
Global Const APPARATUS_SET_B_INPUT = &H9F
Global Const APPARATUS_HOUSE_LIGHT_TIME = 1000
Global Const APPARATUS_HOPPER_TIME = 200
Global Const APPARATUS_OFF_TIME = 100
Global Const APPARATUS_OFF = 0
