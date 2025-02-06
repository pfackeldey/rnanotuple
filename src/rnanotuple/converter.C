#include <vector>
#include <map>

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
        outputFilename = inputFilename;
        outputFilename.ReplaceAll(".root", "_rntuple.root");
    }

    auto inputFile = TFile::Open(inputFilename);

    auto eventsTree = inputFile->Get<TTree>("Events");

    // Loop over all branches and categorize them
    auto branchNames = std::vector<TString>();
    auto countBranchNames = std::vector<TString>();
    auto collectionFields = std::map<TString, std::vector<TString>>();
    auto collectionMaxSizes = std::map<TString, std::size_t>();
    auto recordFields = std::map<TString, std::vector<TString>>();
    auto independentFieldNames = std::vector<TString>();
    for (auto branchPtr : *eventsTree->GetListOfBranches()) {
        auto branch = static_cast<TBranch *>(branchPtr);
        auto branchName = TString(branch->GetName());
        branchNames.push_back(branchName);

        auto firstLeaf = static_cast<TLeaf *>(branch->GetListOfLeaves()->First());
        if (firstLeaf->IsRange()) {
            R__ASSERT(branchName.BeginsWith("n")); // in nanoAOD all count branches start with "n"
            countBranchNames.push_back(branchName);
            auto collectionName = branchName(1, branchName.Length());
            collectionFields[collectionName] = std::vector<TString>();
            collectionMaxSizes[collectionName] = firstLeaf->GetMaximum();
            continue;
        }

        auto leafCount = firstLeaf->GetLeafCount();
        if (leafCount) {
            auto collectionName = TString(leafCount->GetName());
            collectionName = collectionName(1, collectionName.Length());
            R__ASSERT(branchName.BeginsWith(collectionName + "_"));
            R__ASSERT(collectionFields.find(collectionName) != collectionFields.end());
            auto subfieldName = branchName(collectionName.Length()+1, branchName.Length());
            collectionFields[collectionName].push_back(subfieldName);
            continue;
        }

        if (branchName.Contains("_")) {
            auto recordName = branchName(0, branchName.First("_"));
            if (recordFields.find(recordName) == recordFields.end()) {
                recordFields[recordName] = std::vector<TString>();
            }
            recordFields[recordName].push_back(branchName(recordName.Length()+1, branchName.Length()));
            continue;
        }

        independentFieldNames.push_back(branchName);
    }

    if (verbose) {
        // Print the names of attributes, groups, and independent branches
        std::cout << "Collections:" << std::endl;
        for (auto const& [collectionName, subfieldNames] : collectionFields) {
            std::cout << "  " << collectionName << ":" << std::endl;
            for (auto subfieldName : subfieldNames) {
                std::cout << "    " << subfieldName << std::endl;
            }
        }
        std::cout << "Record fields:" << std::endl;
        for (auto const& [recordName, subfieldNames] : recordFields) {
            std::cout << "  " << recordName << ":" << std::endl;
            for (auto subfieldName : subfieldNames) {
                std::cout << "    " << subfieldName << std::endl;
            }
        }
        std::cout << "Independent fields:" << std::endl;
        for (auto independentFieldName : independentFieldNames) {
            std::cout << "    " << independentFieldName << std::endl;
        }
    }

    // Start constructing the RNTuple model
    auto model = RNTupleModel::Create();

    // Independent fields
    for (auto independentFieldName : independentFieldNames) {
        auto branch = eventsTree->GetBranch(independentFieldName);
        auto leaf = branch->GetLeaf(independentFieldName);
        auto type = leaf->GetTypeName();
        auto field = RFieldBase::Create(independentFieldName.Data(), type).Unwrap();
        model->AddField(std::move(field));
    }

    // Record fields
    auto recordOffsets = std::map<TString, std::vector<std::size_t>>();
    for (auto const& [recordName, subfieldNames] : recordFields) {
        std::vector<std::unique_ptr<RFieldBase>> subfields;
        for (auto subfieldName : subfieldNames) {
            auto branch = eventsTree->GetBranch(recordName + "_" + subfieldName);
            auto leaf = branch->GetLeaf(recordName + "_" + subfieldName);
            auto type = leaf->GetTypeName();
            auto field = RFieldBase::Create(subfieldName.Data(), type).Unwrap();
            subfields.push_back(std::move(field));
        }
        auto recordField = std::make_unique<RRecordField>(recordName.Data(), std::move(subfields));
        recordOffsets[recordName] = recordField->GetOffsets();
        model->AddField(std::move(recordField));
    }

    // Collections
    auto collectionOffsets = std::map<TString, std::vector<std::size_t>>();
    auto collectionSizes = std::map<TString, std::size_t>();
    auto collectionLeafSizes = std::map<TString, std::vector<std::size_t>>();
    auto collectionBranchBuffers = std::map<TString, std::vector<std::unique_ptr<unsigned char[]>>>();
    for (auto const& [collectionName, subfieldNames] : collectionFields) {
        std::vector<std::unique_ptr<RFieldBase>> subfields;
        std::vector<TBranch*> subfieldBranches;
        collectionBranchBuffers[collectionName] = std::vector<std::unique_ptr<unsigned char[]>>();
        for (auto subfieldName : subfieldNames) {
            std::cout << collectionName + "_" + subfieldName << std::endl;
            auto branch = eventsTree->GetBranch(collectionName + "_" + subfieldName);
            auto leaf = branch->GetLeaf(collectionName + "_" + subfieldName);
            auto type = leaf->GetTypeName();
            auto field = RFieldBase::Create(subfieldName.Data(), type).Unwrap();
            auto branchBufferSize = collectionMaxSizes[collectionName] * field->GetValueSize();
            subfields.push_back(std::move(field));
            subfieldBranches.push_back(branch);
            auto branchBuffer = std::make_unique<unsigned char[]>(branchBufferSize);
            eventsTree->SetBranchAddress(branch->GetName(), reinterpret_cast<void *>(branchBuffer.get()));
            collectionBranchBuffers[collectionName].emplace_back(std::move(branchBuffer));
        }
        auto recordField = std::make_unique<RRecordField>("_0", std::move(subfields));
        collectionOffsets[collectionName] = recordField->GetOffsets();
        collectionSizes[collectionName] = recordField->GetValueSize();
        auto leafSizes = std::vector<std::size_t>();
        for (std::size_t i = 0; i < recordField->GetSubFields().size(); i++) {
            leafSizes.push_back(recordField->GetSubFields()[i]->GetValueSize());
        }
        collectionLeafSizes[collectionName] = std::move(leafSizes);
        auto collectionField = RVectorField::CreateUntyped(collectionName.Data(), std::move(recordField));
        model->AddField(std::move(collectionField));
    }

    //model->Freeze();
    //auto entry = model->CreateBareEntry();

    // Create the RNTuple writer
    auto writer = RNTupleWriter::Recreate(std::move(model), "Events", outputFilename.Data());

    auto entry = writer->CreateEntry();


    // Bind pointers to the independent fields
    for (auto independentFieldName : independentFieldNames) {
        auto branch = eventsTree->GetBranch(independentFieldName);
        auto leaf = branch->GetLeaf(independentFieldName);
        entry->BindRawPtr(independentFieldName.Data(), leaf->GetValuePointer());
    }

    // Bind pointers to the record subfields
    for (auto const& [recordName, subfieldNames] : recordFields) {
        for (size_t j = 0; j < subfieldNames.size(); j++) {
            auto subfieldBranchName = recordName + "_" + subfieldNames[j];
            auto branch = eventsTree->GetBranch(subfieldBranchName);
            auto leaf = branch->GetLeaf(subfieldBranchName);
            void* fieldPtr = static_cast<char*>(entry->GetPtr<void>(recordName).get()) + recordOffsets[recordName][j];
            eventsTree->SetBranchAddress(subfieldBranchName, fieldPtr);
        }
    }

    // Collections (these are still pretty clunky)
    struct CollectionInfo {
        Int_t maxLength = 0;
        std::unique_ptr<Int_t> countVal;
        std::vector<unsigned char> fieldBuffer;
        std::size_t nLeafs;
        std::size_t recordSize;
        std::vector<std::size_t> offsets;
        std::vector<std::size_t> leafSizes;
        std::vector<std::unique_ptr<unsigned char[]>> leafBuffers;
    };
    auto leafCountCollections = std::map<TString, CollectionInfo>();
    for (auto const& [collectionName, subfieldNames] : collectionFields) {
        auto countBranch = eventsTree->GetBranch("n" + collectionName);
        auto firstLeaf = static_cast<TLeaf *>(countBranch->GetListOfLeaves()->First());
        CollectionInfo c;
        c.maxLength = firstLeaf->GetMaximum();
        c.countVal = std::make_unique<Int_t>();
        c.fieldBuffer.reserve(c.maxLength * collectionSizes[collectionName]);
        c.nLeafs = subfieldNames.size();
        c.recordSize = collectionSizes[collectionName];
        c.offsets = collectionOffsets[collectionName];
        c.leafSizes = collectionLeafSizes[collectionName];
        c.leafBuffers = std::move(collectionBranchBuffers[collectionName]);
        eventsTree->SetBranchAddress(countBranch->GetName(), static_cast<void *>(c.countVal.get()));
        entry->BindRawPtr<void>(collectionName.Data(), &c.fieldBuffer);
        leafCountCollections.emplace(collectionName, std::move(c));
    }


    // Loop over the entries and fill the RNTuple
    auto nEntries = eventsTree->GetEntries();
    for (auto iEntry = 0; iEntry < nEntries; iEntry++) {
        if (verbose && iEntry % 1000 == 0) {
            std::cout << "Processing entry " << iEntry << " of " << nEntries << std::endl;
        }
        eventsTree->GetEntry(iEntry);

        for (auto &[collectionName, c] : leafCountCollections) {
            const auto sizeOfRecord = c.recordSize;
            std::cout << collectionName << " Resizing to " << sizeOfRecord * (*c.countVal) << std::endl;
            c.fieldBuffer.resize(sizeOfRecord * (*c.countVal));

            const auto nLeafs = c.nLeafs;
            for (std::size_t l = 0; l < nLeafs; ++l) {
                const auto offset = c.offsets[l];
                const auto sizeOfLeaf = c.leafSizes[l];
                const auto &leafBuffer = c.leafBuffers[l];
                for (Int_t j = 0; j < *c.countVal; ++j) {
                    std::cout << "Copying " << j * sizeOfRecord + offset << " from " << j * sizeOfLeaf << " size " << sizeOfLeaf << std::endl;
                    std::memcpy(c.fieldBuffer.data() + j * sizeOfRecord + offset, leafBuffer.get() + (j * sizeOfLeaf), sizeOfLeaf);
                    std::cout << "Source data:" << *(int8_t*)(leafBuffer.get() + (j * sizeOfLeaf)) << std::endl;
                    std::cout << "Dest data:" << *(int8_t*)(c.fieldBuffer.data() + j * sizeOfRecord + offset) << std::endl;
                }
            }
        }

        writer->Fill(*entry);
    }

}
