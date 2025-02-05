#include <vector>
#include <set>
#include <cassert>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TLeaf.h"
#include <ROOT/RNTuple.hxx>
#include <ROOT/RFieldBase.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RField/RFieldRecord.hxx>
#include <ROOT/RField/RFieldSequenceContainer.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>


using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;
using RFieldBase = ROOT::Experimental::RFieldBase;
using RRecordField = ROOT::Experimental::RRecordField;
using RVectorField = ROOT::Experimental::RVectorField;


void converter(TString inputFilename, TString outputFilename = "", Bool_t verbose = kTRUE) {
    if (outputFilename == "") {
        outputFilename = inputFilename(0, inputFilename.Length() - 5) + "_rntuple.root";
    }

    auto inputFile = TFile::Open(inputFilename);

    auto eventsTree = inputFile->Get<TTree>("Events");

    // Names of all branches in the tree
    auto branchNames = std::vector<TString>();
    for (auto branch : *eventsTree->GetListOfBranches()) {
        branchNames.push_back(branch->GetName());
    }

    // Names of branches describing the jaggedness of arrays
    // They are assumed to be all that start with "n"
    auto jaggednessBranchNames = std::vector<TString>();
    for (auto branchName : branchNames) {
        if (branchName.BeginsWith("n")) {
            jaggednessBranchNames.push_back(branchName(1, branchName.Length() - 1));
        }
    }

    // Names of branches associated with jaggedness branches
    auto attributeBranchNames = std::vector<std::vector<TString>>(jaggednessBranchNames.size());
    for (auto branchName : branchNames) {
        for (size_t i = 0; i < jaggednessBranchNames.size(); i++) {
            if (branchName.BeginsWith(jaggednessBranchNames[i])) {
                attributeBranchNames[i].push_back(branchName);
            }
        }
    }
    auto allAttributeBranchNames = std::vector<TString>();
    for (auto attributeBranchNames_ : attributeBranchNames) {
        allAttributeBranchNames.insert(allAttributeBranchNames.end(), attributeBranchNames_.begin(), attributeBranchNames_.end());
    }

    // Names of groups of branches that are not jagged
    auto groupNamesSet = std::set<TString>();
    for (auto branchName : branchNames) {
        if (branchName.BeginsWith("n") || !branchName.Contains("_")) {
            continue;
        }
        // Make sure they are not already in attributeBranchNames
        if (std::find(allAttributeBranchNames.begin(), allAttributeBranchNames.end(), branchName) != allAttributeBranchNames.end()) {
            continue;
        }
        groupNamesSet.insert(branchName(0, branchName.First("_")));
    }
    auto groupNames = std::vector<TString>(groupNamesSet.begin(), groupNamesSet.end());

    // Names of the branches in each group
    auto groupBranchNames = std::vector<std::vector<TString>>(groupNames.size());
    for (auto branchName : branchNames) {
        for (size_t i = 0; i < groupNames.size(); i++) {
            if (branchName.BeginsWith(groupNames[i]+"_")) {
                groupBranchNames[i].push_back(branchName);
            }
        }
    }
    auto allGroupBranchNames = std::vector<TString>();
    for (auto groupBranchNames_ : groupBranchNames) {
        allGroupBranchNames.insert(allGroupBranchNames.end(), groupBranchNames_.begin(), groupBranchNames_.end());
    }

    // Independent branches
    auto independentBranchNames = std::vector<TString>();
    for (auto branchName : branchNames) {
        if (branchName.BeginsWith("n")) {
            continue;
        }
        if (std::find(allAttributeBranchNames.begin(), allAttributeBranchNames.end(), branchName) != allAttributeBranchNames.end()) {
            continue;
        }
        if (std::find(allGroupBranchNames.begin(), allGroupBranchNames.end(), branchName) != allGroupBranchNames.end()) {
            continue;
        }
        independentBranchNames.push_back(branchName);
    }

    if (verbose) {
        // Print the names of attributes, groups, and independent branches
        std::cout << "Attributes:" << std::endl;
        for (size_t i = 0; i < jaggednessBranchNames.size(); i++) {
            std::cout << jaggednessBranchNames[i] << ":" << std::endl;
            for (auto attributeBranchName : attributeBranchNames[i]) {
                std::cout << "    " << attributeBranchName << std::endl;
            }
        }
        std::cout << "Groups:" << std::endl;
        for (size_t i = 0; i < groupNames.size(); i++) {
            std::cout << groupNames[i] << ":" << std::endl;
            for (auto groupBranchName : groupBranchNames[i]) {
                std::cout << "    " << groupBranchName << std::endl;
            }
        }
        std::cout << "Independent:" << std::endl;
        for (auto independentBranchName : independentBranchNames) {
            std::cout << "    " << independentBranchName << std::endl;
        }
    }

    // Start constructing the RNTuple model
    auto model = RNTupleModel::Create();

    // Independent branches
    for (auto independentBranchName : independentBranchNames) {
        auto branch = eventsTree->GetBranch(independentBranchName);
        auto leaf = branch->GetLeaf(independentBranchName);
        auto type = leaf->GetTypeName();
        auto field = RFieldBase::Create(independentBranchName.Data(), type).Unwrap();
        model->AddField(std::move(field));
    }

    // Jagged branches
    for (size_t i = 0; i < jaggednessBranchNames.size(); i++) {
        auto jaggednessBranchName = jaggednessBranchNames[i];
        std::vector<std::unique_ptr<RFieldBase>> leafFields;
        for (auto attributeBranchName : attributeBranchNames[i]) {
            auto branch = eventsTree->GetBranch(attributeBranchName);
            auto leaf = branch->GetLeaf(attributeBranchName);
            auto type = leaf->GetTypeName();
            auto field = RFieldBase::Create(attributeBranchName(jaggednessBranchName.Length()+1, attributeBranchName.Length()-1), type).Unwrap();
            leafFields.push_back(std::move(field));
        }
        auto recordField = std::make_unique<RRecordField>("_0", std::move(leafFields));
        auto recordFieldPtr = recordField.get();
        auto collectionField = RVectorField::CreateUntyped(jaggednessBranchName.Data(), std::move(recordField));
        model->AddField(std::move(collectionField));
        // Add projected subfields for leafs
        for (auto leaf : recordFieldPtr->GetSubFields()) {
            auto name = leaf->GetFieldName();
            auto fullname = std::string(jaggednessBranchName.Data()) + "_" + name;
            //std::cout << "fullname: " << fullname << std::endl;
            auto type = leaf->GetTypeName();
            //std::cout << name << " " << jaggednessBranchName << std::endl;
            //std::cout << "Search: " << model->FindField(jaggednessBranchName) << std::endl;
            auto projectedField = RFieldBase::Create(fullname, "ROOT::VecOps::RVec<" + type + ">").Unwrap();
            //std::cout << "projectedField: " << projectedField->GetFieldName() << std::endl;
            std::string jaggednessBranchName_ = jaggednessBranchName.Data();
            model->AddProjectedField(std::move(projectedField), [&jaggednessBranchName_, &name, &fullname](const std::string &fieldName) {
                //cout << fieldName << " " << jaggednessBranchName_ << " " << jaggednessBranchName_ + "._0." + name << endl;
                if (fieldName == fullname)
                    return jaggednessBranchName_;
                else
                    return jaggednessBranchName_ + "._0." + name;
            });
        }

    }

    // Groups
    std::vector<std::vector<std::size_t>> groupOffsets;
    std::vector<std::size_t> groupSizes;
    for (size_t i = 0; i < groupNames.size(); i++) {
        auto groupName = groupNames[i];
        std::vector<std::unique_ptr<RFieldBase>> leafFields;
        for (auto groupBranchName : groupBranchNames[i]) {
            auto branch = eventsTree->GetBranch(groupBranchName);
            auto leaf = branch->GetLeaf(groupBranchName);
            auto type = leaf->GetTypeName();
            auto field = RFieldBase::Create(groupBranchName(groupName.Length()+1, groupBranchName.Length()-1), type).Unwrap();
            leafFields.push_back(std::move(field));
        }
        auto recordField = std::make_unique<RRecordField>(groupName.Data(), std::move(leafFields));
        groupOffsets.push_back(recordField->GetOffsets());
        groupSizes.push_back(recordField->GetValueSize());
        model->AddField(std::move(recordField));
    }

    auto& defEntry = model->GetDefaultEntry();

    // Bind pointers to the independent leaves
    for (auto independentBranchName : independentBranchNames) {
        auto branch = eventsTree->GetBranch(independentBranchName);
        auto leaf = branch->GetLeaf(independentBranchName);
        defEntry.BindRawPtr(independentBranchName.Data(), leaf->GetValuePointer());
    }

    // Jagged branches
    for (size_t i = 0; i < jaggednessBranchNames.size(); i++) {
        auto jaggednessBranchName = jaggednessBranchNames[i];
        auto nBranch = eventsTree->GetBranch("n" + jaggednessBranchName);
        auto nLeaf = nBranch->GetLeaf("n" + jaggednessBranchName);
        //auto n = nLeaf->GetValue<int>();
        for (size_t j = 0; j<attributeBranchNames[i].size(); j++) {
            auto attributeBranchName = attributeBranchNames[i][j];
            auto branch = eventsTree->GetBranch(attributeBranchName);
            auto leaf = branch->GetLeaf(attributeBranchName);
            //std::cout << jaggednessBranchName << " " << attributeBranchName << " " << leaf->GetTypeName() << std::endl;
            //defEntry.addValue();
            //std::cout << defEntry.GetPtr<void>(jaggednessBranchName + "._0." + attributeBranchName(jaggednessBranchName.Length()+1, attributeBranchName.Length())).get() << std::endl;
            //void* fieldPtr = static_cast<char*>(defEntry.GetPtr<void>(groupName).get()) + groupOffsets[i][j];
            //std::cout << fieldPtr << std::endl;
            //eventsTree->SetBranchAddress(groupBranchName, fieldPtr);
            //break;
            //entry->BindRawPtr((groupName + "." + groupBranchName(groupName.Length()+1, groupBranchName.Length()-groupName.Length()-1)).Data(), leaf->GetValuePointer());
        }
    }

    // Bind pointers to the group leaves
    for (size_t i = 0; i < groupNames.size(); i++) {
        auto groupName = groupNames[i];
        for (size_t j = 0; j<groupBranchNames[i].size(); j++) {
            auto groupBranchName = groupBranchNames[i][j];
            auto branch = eventsTree->GetBranch(groupBranchName);
            auto leaf = branch->GetLeaf(groupBranchName);
            //std::cout << groupName << " " << groupBranchName << " " << leaf->GetTypeName() << " " << groupSizes[i] << std::endl;
            //std::cout << defEntry.GetPtr<void>(groupName).get() << std::endl;
            void* fieldPtr = static_cast<char*>(defEntry.GetPtr<void>(groupName).get()) + groupOffsets[i][j];
            //std::cout << fieldPtr << std::endl;
            eventsTree->SetBranchAddress(groupBranchName, fieldPtr);
            //break;
            //entry->BindRawPtr((groupName + "." + groupBranchName(groupName.Length()+1, groupBranchName.Length()-groupName.Length()-1)).Data(), leaf->GetValuePointer());
        }
    }

    // Create the RNTuple writer
    auto writer = RNTupleWriter::Recreate(std::move(model), "Events", outputFilename.Data());

    auto nEntries = eventsTree->GetEntries();

    //auto entry = writer->CreateEntry();
    //std::cout << entry->GetPtr<void>("Electronss").get() << std::endl;




    // Loop over the entries and fill the RNTuple
    for (auto iEntry = 0; iEntry < nEntries; iEntry++) {
        if (verbose && iEntry % 1000 == 0) {
            std::cout << "Processing entry " << iEntry << " of " << nEntries << std::endl;
        }
        eventsTree->GetEntry(iEntry);

        writer->Fill();
    }

}
