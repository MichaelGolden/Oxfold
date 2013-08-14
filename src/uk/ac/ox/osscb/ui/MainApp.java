package uk.ac.ox.osscb.ui;

import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.swing.JLabel;
import javax.swing.SwingConstants;
import java.awt.FlowLayout;
import javax.swing.JButton;
import java.awt.BorderLayout;
import java.awt.GridBagLayout;
import javax.swing.JRadioButton;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import java.awt.Component;
import javax.swing.JCheckBox;
import javax.swing.JTextField;
import java.awt.GridLayout;

public class MainApp {

	private JFrame frmOxfold;
	private JTextField textField;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					MainApp window = new MainApp();
					window.frmOxfold.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public MainApp() {
		initialize();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {
		frmOxfold = new JFrame();
		frmOxfold.setTitle("Oxfold");
		frmOxfold.setBounds(100, 100, 626, 416);
		frmOxfold.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frmOxfold.getContentPane().setLayout(new BoxLayout(frmOxfold.getContentPane(), BoxLayout.X_AXIS));
		
		JPanel panel_4 = new JPanel();
		frmOxfold.getContentPane().add(panel_4);
		panel_4.setLayout(new BoxLayout(panel_4, BoxLayout.Y_AXIS));
		
		JPanel panel_1 = new JPanel();
		panel_4.add(panel_1);
		panel_1.setBorder(new TitledBorder(null, "Select an alignment file", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		panel_1.setLayout(new BorderLayout(0, 0));
		
		JPanel panel_3 = new JPanel();
		panel_1.add(panel_3, BorderLayout.CENTER);
		GridBagLayout gbl_panel_3 = new GridBagLayout();
		gbl_panel_3.columnWidths = new int[]{0, 0, 0, 0};
		gbl_panel_3.rowHeights = new int[]{0, 0};
		gbl_panel_3.columnWeights = new double[]{0.0, 0.0, 0.0, Double.MIN_VALUE};
		gbl_panel_3.rowWeights = new double[]{0.0, Double.MIN_VALUE};
		panel_3.setLayout(gbl_panel_3);
		
		JLabel lblDropAndDrag = new JLabel("Drop and drag a file or ");
		GridBagConstraints gbc_lblDropAndDrag = new GridBagConstraints();
		gbc_lblDropAndDrag.insets = new Insets(0, 0, 0, 5);
		gbc_lblDropAndDrag.gridx = 0;
		gbc_lblDropAndDrag.gridy = 0;
		panel_3.add(lblDropAndDrag, gbc_lblDropAndDrag);
		
		JButton btnNewButton = new JButton("Browse..");
		GridBagConstraints gbc_btnNewButton = new GridBagConstraints();
		gbc_btnNewButton.insets = new Insets(0, 0, 0, 5);
		gbc_btnNewButton.gridx = 1;
		gbc_btnNewButton.gridy = 0;
		panel_3.add(btnNewButton, gbc_btnNewButton);
		
		JLabel lblmostFormatsSupported = new JLabel("(most formats supported)");
		GridBagConstraints gbc_lblmostFormatsSupported = new GridBagConstraints();
		gbc_lblmostFormatsSupported.gridx = 2;
		gbc_lblmostFormatsSupported.gridy = 0;
		panel_3.add(lblmostFormatsSupported, gbc_lblmostFormatsSupported);
		
		JPanel panel = new JPanel();
		panel_4.add(panel);
		panel.setBorder(new TitledBorder(null, "Select a phylogenetic tree", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagLayout gbl_panel = new GridBagLayout();
		gbl_panel.columnWidths = new int[]{0, 0, 0};
		gbl_panel.rowHeights = new int[]{0, 0, 0, 0};
		gbl_panel.columnWeights = new double[]{0.0, 0.0, Double.MIN_VALUE};
		gbl_panel.rowWeights = new double[]{0.0, 0.0, 0.0, Double.MIN_VALUE};
		panel.setLayout(gbl_panel);
		
		JRadioButton rdbtnInfer = new JRadioButton("Infer tree using FastTree (recommended)");
		GridBagConstraints gbc_rdbtnInfer = new GridBagConstraints();
		gbc_rdbtnInfer.insets = new Insets(0, 0, 5, 5);
		gbc_rdbtnInfer.anchor = GridBagConstraints.WEST;
		gbc_rdbtnInfer.gridx = 0;
		gbc_rdbtnInfer.gridy = 0;
		panel.add(rdbtnInfer, gbc_rdbtnInfer);
		
		JRadioButton rdbtnNewRadioButton = new JRadioButton("Specify a newick tree file");
		GridBagConstraints gbc_rdbtnNewRadioButton = new GridBagConstraints();
		gbc_rdbtnNewRadioButton.insets = new Insets(0, 0, 5, 5);
		gbc_rdbtnNewRadioButton.anchor = GridBagConstraints.WEST;
		gbc_rdbtnNewRadioButton.gridx = 0;
		gbc_rdbtnNewRadioButton.gridy = 1;
		panel.add(rdbtnNewRadioButton, gbc_rdbtnNewRadioButton);
		
		JButton btnBrowse_1 = new JButton("Browse..");
		GridBagConstraints gbc_btnBrowse_1 = new GridBagConstraints();
		gbc_btnBrowse_1.insets = new Insets(0, 0, 5, 0);
		gbc_btnBrowse_1.gridx = 1;
		gbc_btnBrowse_1.gridy = 1;
		panel.add(btnBrowse_1, gbc_btnBrowse_1);
		
		JRadioButton rdbtnUseNonevolutionaryModel = new JRadioButton("Use non-evolutionary model (no tree)");
		GridBagConstraints gbc_rdbtnUseNonevolutionaryModel = new GridBagConstraints();
		gbc_rdbtnUseNonevolutionaryModel.anchor = GridBagConstraints.WEST;
		gbc_rdbtnUseNonevolutionaryModel.insets = new Insets(0, 0, 0, 5);
		gbc_rdbtnUseNonevolutionaryModel.gridx = 0;
		gbc_rdbtnUseNonevolutionaryModel.gridy = 2;
		panel.add(rdbtnUseNonevolutionaryModel, gbc_rdbtnUseNonevolutionaryModel);
		
		JPanel panel_2 = new JPanel();
		panel_4.add(panel_2);
		panel_2.setBorder(new TitledBorder(null, "Options", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagLayout gbl_panel_2 = new GridBagLayout();
		gbl_panel_2.columnWidths = new int[]{0, 0, 0};
		gbl_panel_2.rowHeights = new int[]{0, 0};
		gbl_panel_2.columnWeights = new double[]{0.0, 0.0, Double.MIN_VALUE};
		gbl_panel_2.rowWeights = new double[]{0.0, Double.MIN_VALUE};
		panel_2.setLayout(gbl_panel_2);
		
		JLabel lblNumberOfProcessors = new JLabel("Number of processors");
		GridBagConstraints gbc_lblNumberOfProcessors = new GridBagConstraints();
		gbc_lblNumberOfProcessors.insets = new Insets(0, 0, 0, 5);
		gbc_lblNumberOfProcessors.gridx = 0;
		gbc_lblNumberOfProcessors.gridy = 0;
		panel_2.add(lblNumberOfProcessors, gbc_lblNumberOfProcessors);
		
		JSpinner spinner = new JSpinner();
		spinner.setModel(new SpinnerNumberModel(new Integer(1), null, null, new Integer(1)));
		GridBagConstraints gbc_spinner = new GridBagConstraints();
		gbc_spinner.gridx = 1;
		gbc_spinner.gridy = 0;
		panel_2.add(spinner, gbc_spinner);
		
		JPanel panel_5 = new JPanel();
		panel_5.setBorder(new TitledBorder(null, "Output", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		frmOxfold.getContentPane().add(panel_5);
		panel_5.setLayout(new BorderLayout(0, 0));
		
		JPanel panel_6 = new JPanel();
		panel_5.add(panel_6, BorderLayout.CENTER);
		panel_6.setLayout(new BoxLayout(panel_6, BoxLayout.Y_AXIS));
		
		JPanel panel_7 = new JPanel();
		panel_6.add(panel_7);
		GridBagLayout gbl_panel_7 = new GridBagLayout();
		gbl_panel_7.columnWidths = new int[]{0, 0, 0};
		gbl_panel_7.rowHeights = new int[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		gbl_panel_7.columnWeights = new double[]{1.0, 0.0, Double.MIN_VALUE};
		gbl_panel_7.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
		panel_7.setLayout(gbl_panel_7);
		
		JLabel lblNewLabel = new JLabel("Structure outputs");
		GridBagConstraints gbc_lblNewLabel = new GridBagConstraints();
		gbc_lblNewLabel.anchor = GridBagConstraints.WEST;
		gbc_lblNewLabel.insets = new Insets(0, 0, 5, 5);
		gbc_lblNewLabel.gridx = 0;
		gbc_lblNewLabel.gridy = 0;
		panel_7.add(lblNewLabel, gbc_lblNewLabel);
		
		JCheckBox chckbxConnectFormatct = new JCheckBox("Connect format (.ct)");
		GridBagConstraints gbc_chckbxConnectFormatct = new GridBagConstraints();
		gbc_chckbxConnectFormatct.anchor = GridBagConstraints.WEST;
		gbc_chckbxConnectFormatct.insets = new Insets(0, 0, 5, 5);
		gbc_chckbxConnectFormatct.gridx = 0;
		gbc_chckbxConnectFormatct.gridy = 1;
		panel_7.add(chckbxConnectFormatct, gbc_chckbxConnectFormatct);
		
		JCheckBox chckbxViennaDotBracket = new JCheckBox("Vienna dot bracket (.dbn)");
		GridBagConstraints gbc_chckbxViennaDotBracket = new GridBagConstraints();
		gbc_chckbxViennaDotBracket.anchor = GridBagConstraints.WEST;
		gbc_chckbxViennaDotBracket.insets = new Insets(0, 0, 5, 5);
		gbc_chckbxViennaDotBracket.gridx = 0;
		gbc_chckbxViennaDotBracket.gridy = 2;
		panel_7.add(chckbxViennaDotBracket, gbc_chckbxViennaDotBracket);
		
		JLabel lblAdditionalOutputs = new JLabel("Additional outputs");
		GridBagConstraints gbc_lblAdditionalOutputs = new GridBagConstraints();
		gbc_lblAdditionalOutputs.anchor = GridBagConstraints.WEST;
		gbc_lblAdditionalOutputs.insets = new Insets(0, 0, 5, 5);
		gbc_lblAdditionalOutputs.gridx = 0;
		gbc_lblAdditionalOutputs.gridy = 3;
		panel_7.add(lblAdditionalOutputs, gbc_lblAdditionalOutputs);
		
		JCheckBox chckbxSummaryFile = new JCheckBox("Summary file (.txt)");
		GridBagConstraints gbc_chckbxSummaryFile = new GridBagConstraints();
		gbc_chckbxSummaryFile.anchor = GridBagConstraints.WEST;
		gbc_chckbxSummaryFile.insets = new Insets(0, 0, 5, 5);
		gbc_chckbxSummaryFile.gridx = 0;
		gbc_chckbxSummaryFile.gridy = 4;
		panel_7.add(chckbxSummaryFile, gbc_chckbxSummaryFile);
		
		JCheckBox chckbxBasepairProbabilityMatrices = new JCheckBox("Base-pair probability matrix (.bp)");
		GridBagConstraints gbc_chckbxBasepairProbabilityMatrices = new GridBagConstraints();
		gbc_chckbxBasepairProbabilityMatrices.anchor = GridBagConstraints.WEST;
		gbc_chckbxBasepairProbabilityMatrices.insets = new Insets(0, 0, 5, 5);
		gbc_chckbxBasepairProbabilityMatrices.gridx = 0;
		gbc_chckbxBasepairProbabilityMatrices.gridy = 5;
		panel_7.add(chckbxBasepairProbabilityMatrices, gbc_chckbxBasepairProbabilityMatrices);
		
		JCheckBox chckbxCircularPlotpng = new JCheckBox("Circular plot (circular.png)");
		GridBagConstraints gbc_chckbxCircularPlotpng = new GridBagConstraints();
		gbc_chckbxCircularPlotpng.anchor = GridBagConstraints.WEST;
		gbc_chckbxCircularPlotpng.insets = new Insets(0, 0, 5, 5);
		gbc_chckbxCircularPlotpng.gridx = 0;
		gbc_chckbxCircularPlotpng.gridy = 6;
		panel_7.add(chckbxCircularPlotpng, gbc_chckbxCircularPlotpng);
		
		JCheckBox chckbxCircularPlotsvg = new JCheckBox("Circular plot (circular.svg)");
		GridBagConstraints gbc_chckbxCircularPlotsvg = new GridBagConstraints();
		gbc_chckbxCircularPlotsvg.anchor = GridBagConstraints.WEST;
		gbc_chckbxCircularPlotsvg.insets = new Insets(0, 0, 5, 5);
		gbc_chckbxCircularPlotsvg.gridx = 0;
		gbc_chckbxCircularPlotsvg.gridy = 7;
		panel_7.add(chckbxCircularPlotsvg, gbc_chckbxCircularPlotsvg);
		
		JCheckBox chckbxStructurePlotpng = new JCheckBox("Structure plot (_structure.png)");
		GridBagConstraints gbc_chckbxStructurePlotpng = new GridBagConstraints();
		gbc_chckbxStructurePlotpng.insets = new Insets(0, 0, 5, 5);
		gbc_chckbxStructurePlotpng.anchor = GridBagConstraints.WEST;
		gbc_chckbxStructurePlotpng.gridx = 0;
		gbc_chckbxStructurePlotpng.gridy = 8;
		panel_7.add(chckbxStructurePlotpng, gbc_chckbxStructurePlotpng);
		
		JLabel lblExportPath = new JLabel("Export path");
		GridBagConstraints gbc_lblExportPath = new GridBagConstraints();
		gbc_lblExportPath.anchor = GridBagConstraints.WEST;
		gbc_lblExportPath.insets = new Insets(0, 0, 5, 5);
		gbc_lblExportPath.gridx = 0;
		gbc_lblExportPath.gridy = 10;
		panel_7.add(lblExportPath, gbc_lblExportPath);
		
		textField = new JTextField();
		GridBagConstraints gbc_textField = new GridBagConstraints();
		gbc_textField.gridwidth = 2;
		gbc_textField.insets = new Insets(0, 0, 5, 0);
		gbc_textField.fill = GridBagConstraints.HORIZONTAL;
		gbc_textField.gridx = 0;
		gbc_textField.gridy = 11;
		panel_7.add(textField, gbc_textField);
		textField.setColumns(10);
		
		JButton btnNewButton_1 = new JButton("Export");
		GridBagConstraints gbc_btnNewButton_1 = new GridBagConstraints();
		gbc_btnNewButton_1.gridx = 1;
		gbc_btnNewButton_1.gridy = 12;
		panel_7.add(btnNewButton_1, gbc_btnNewButton_1);
	}

}
