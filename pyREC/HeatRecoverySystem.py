from AuxiliaryComponent import AuxiliaryComponent


class HeatRecoverySystem(AuxiliaryComponent):
    def __init__(self, id, eff, replacement_year, replacement_cost, tech='HeatRecoverySystem', ):
        """

        :param id:
        :param eff:
        :param replacement_year:
        :param replacement_cost: initial cost [€]
        :param tech:

        """
        super().__init__(id=id, tech=tech, replacement_cost=replacement_cost, replacement_year=replacement_year)

        self.eff = eff

    def compute_output(self, power_heat_input):
        """

        :param power_heat_input: [kW]
        :return: power_heat_output [kW]
        """
        power_heat_output = self.eff * power_heat_input

        return power_heat_input
